from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.ticker import ScalarFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
import argparse

def parse_file(path):
    meta = {} # metadata
    snapshots = [] # data blocks

    with open(path, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    i = 0 # line counter

    # parse header, extract metadata
    while i < len(lines) and lines[i].startswith("#"):
        line = lines[i]
        if ":" in line:
            key, val = line[1:].split(":", 1)
            meta[key.strip()] = val.strip()
        i += 1

    nx = int(meta["nx"])
    ny = int(meta["ny"])

    # parse timesteps
    while i < len(lines):
        if lines[i].startswith("$ timestep"):
            timestep = int(lines[i].split()[-1]) # extract timestep
            i += 1

            t = float(lines[i].split("=")[-1].strip()) # extract time
            i += 1

            block = []
            for _ in range(nx):
                row = [float(x) for x in lines[i].split()] # collect nx lines
                block.append(row)
                i += 1

            snapshots.append({
                "timestep": timestep,
                "t": t,
                "data": np.array(block)
            })
        else:
            i += 1 # handle empty lines

    meta = {
        "rank": int(meta["rank"]),
        "grid_i": int(meta["grid i"]),
        "grid_j": int(meta["grid j"]),
        "i0": int(meta["i0"]),
        "j0": int(meta["j0"]),
        "nx": nx,
        "ny": ny,
    }

    return meta, snapshots

def parse_files(dir):
    paths = sorted(Path(dir).glob("*rank*.dat"))
    return [parse_file(p) for p in paths]

# extract info from global metadata file
def parse_metadata(path):
    path = Path(path)
    meta = {}

    converters = {
        "nx": int,
        "ny": int,
        "lx": float,
        "ly": float,
        "dx": float,
        "dy": float,
        "dt0": float,
        "rho": float,
        "nu": float,
        "courant": float,
        "p_max": float,
        "t_end": float,
        "dt_out": float,
    }

    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()

            if not line or not line.startswith("#"):
                continue

            content = line[1:].strip()
            if ":" not in content:
                continue

            key, value = content.split(":", 1)
            key = key.strip().lower()
            value = value.strip()

            if key in converters:
                meta[key] = converters[key](value)
            else:
                meta[key] = value

    return meta

def stitch_snapshot(parsed_files, snapshot_index, nx_global, ny_global):

    field = np.zeros((nx_global, ny_global))

    for meta, snapshots in parsed_files:
        block = snapshots[snapshot_index]["data"]
        i0, j0 = meta["i0"], meta["j0"] # global start element
        nx, ny = meta["nx"], meta["ny"]

        # stitch block into global domain
        field[i0:i0+nx, j0:j0+ny] = block

    return field

def construct_history(data):
    nx_global = max(meta["i0"] + meta["nx"] for meta, _ in data)
    ny_global = max(meta["j0"] + meta["ny"] for meta, _ in data)
    num_steps = data[0][1][-1]["timestep"] + 1 # get timestep (0 indexed)

    history = np.zeros(shape=(num_steps, nx_global, ny_global))
    for step in range(num_steps):
        history[step] = stitch_snapshot(data, step, nx_global, ny_global) # stitch data blocks together
    return history

# compute limits
def robust_limits(data, symmetric=False, q=(1.0, 99.0)):
    arr = np.asarray(data, dtype=float)
    lo, hi = np.nanpercentile(arr, q)

    if symmetric:
        m = max(abs(lo), abs(hi))
        return -m, m
    return lo, hi

# normalise the data across time
def make_norm(history, field_type="scalar", robust=True):
    if field_type == "diverging":
        vmin, vmax = robust_limits(history, symmetric=True) if robust else (
            float(np.nanmin(history)), float(np.nanmax(history))
        )
        m = max(abs(vmin), abs(vmax))
        return TwoSlopeNorm(vmin=-m, vcenter=0.0, vmax=m) # type: ignore
    else:
        vmin, vmax = robust_limits(history, symmetric=False) if robust else (
            float(np.nanmin(history)), float(np.nanmax(history))
        )
        return Normalize(vmin=vmin, vmax=vmax) # type: ignore


def build_times_from_metadata(meta, nt):
    dt0 = float(meta.get("dt0", 0.0))
    dt_out = float(meta["dt_out"])
    return dt0 + np.arange(nt) * dt_out


def save_triptych_gif(
    P_history,
    u_history,
    v_history,
    gif_path,
    *,
    meta,
    p_norm,
    u_norm,
    v_norm,
    times=None,
    fps=12,
    interval=80,
    dpi=170,
    interpolation="nearest",
    transpose=True,
):
    P_history = np.asarray(P_history, dtype=float)
    u_history = np.asarray(u_history, dtype=float)
    v_history = np.asarray(v_history, dtype=float)

    nt = P_history.shape[0]
    lx = float(meta["lx"])
    ly = float(meta["ly"])
    extent = [0.0, lx, 0.0, ly]

    plt.rcParams.update({
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "axes.edgecolor": "0.15",
        "axes.labelcolor": "0.1",
        "xtick.color": "0.2",
        "ytick.color": "0.2",
        "font.size": 11,
        "axes.titlesize": 13,
        "axes.labelsize": 11,
    })

    fig, axs = plt.subplots(
        1, 3,
        figsize=(18, 6),
        sharex=True,
        sharey=True,
    )
    fig.tight_layout()

    def frame_data(arr, k):
        return arr[k].T if transpose else arr[k]

    imP = axs[0].imshow(
        frame_data(P_history, 0),
        origin="lower",
        extent=extent,
        aspect="equal",
        interpolation=interpolation,
        cmap="cividis",
        norm=p_norm,
    )

    imu = axs[1].imshow(
        frame_data(u_history, 0),
        origin="lower",
        extent=extent,
        aspect="equal",
        interpolation=interpolation,
        cmap="RdBu_r",
        norm=u_norm,
    )

    imv = axs[2].imshow(
        frame_data(v_history, 0),
        origin="lower",
        extent=extent,
        aspect="equal",
        interpolation=interpolation,
        cmap="RdBu_r",
        norm=v_norm,
    )

    axs[0].set_title("Pressure")
    axs[1].set_title("x-velocity")
    axs[2].set_title("y-velocity")

    for i, ax in enumerate(axs):
        ax.set_xlim(0.0, lx)
        ax.set_ylim(0.0, ly)
        ax.set_xlabel("x [m]")
        if i == 0:
            ax.set_ylabel("y [m]")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    def add_horizontal_cbar(fig, ax, im, label):
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("bottom", size="5%", pad=0.75)
        cbar = fig.colorbar(im, cax=cax, orientation="horizontal")
        cbar.set_label(label)
        cbar.ax.xaxis.set_ticks_position("bottom")
        cbar.ax.xaxis.set_label_position("bottom")
        return cbar

    cbarP = add_horizontal_cbar(fig, axs[0], imP, "Pressure [Pa]")
    cbaru = add_horizontal_cbar(fig, axs[1], imu, "u [m s$^{-1}$]")
    cbarv = add_horizontal_cbar(fig, axs[2], imv, "v [m s$^{-1}$]")

    for cbar in (cbarP, cbaru, cbarv):
        cbar.minorticks_on()
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_powerlimits((0, 0))
        cbar.formatter = formatter
        cbar.update_ticks()

    if times is None:
        suptitle = fig.suptitle("Flow fields | step 0", fontsize=16)
    else:
        suptitle = fig.suptitle(f"Flow fields | t = {times[0]:.3f} s", fontsize=16)

    footer = fig.text(
        0.5, 0.02,
        f"Domain: {lx:.3g} m × {ly:.3g} m",
        ha="center", va="bottom",
        fontsize=12, color="0.35"
    )

    def update(frame):
        imP.set_data(frame_data(P_history, frame))
        imu.set_data(frame_data(u_history, frame))
        imv.set_data(frame_data(v_history, frame))

        if times is None:
            suptitle.set_text(f"Flow fields | step {frame}")
        else:
            suptitle.set_text(f"Flow fields | t = {times[frame]:.3f} s")

        return imP, imu, imv, suptitle, footer

    anim = FuncAnimation(
        fig,
        update,
        frames=nt,
        interval=interval,
        blit=False,
        repeat=True,
    )

    gif_path = Path(gif_path)
    gif_path.parent.mkdir(parents=True, exist_ok=True)
    anim.save(gif_path, writer=PillowWriter(fps=fps), dpi=dpi) # type: ignore
    plt.close(fig)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Postprocess CW-MPI output into gifs and diagnostics")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--core", action="store_true", help="Process core output under ./out")
    group.add_argument("--extension", action="store_true", help="Process extension output under ./out/extension")
    parser.add_argument("--root", type=str, help="Path to output root (overrides --core/--extension)")
    args = parser.parse_args()

    if args.root:
        root = Path(args.root)
    elif args.core:
        root = Path("./out/core")
    elif args.extension:
        root = Path("./out/extension")
    else:
        root = Path("./out")

    print(f"Using root: {root}")

    print("Parsing metadata...")
    meta = parse_metadata(root / "metadata.dat")
    print("Done.\n")

    print("Parsing files...")
    P_data = parse_files(root / "P")
    u_data = parse_files(root / "u")
    v_data = parse_files(root / "v")
    print("Done.\n")

    print("Stiching history...")
    P_history = construct_history(P_data)
    u_history = construct_history(u_data)
    v_history = construct_history(v_data)
    print("Done.\n")

    times = build_times_from_metadata(meta, P_history.shape[0])
    p_norm = make_norm(P_history, field_type="scalar", robust=True)
    u_norm = make_norm(u_history, field_type="diverging", robust=True)
    v_norm = make_norm(v_history, field_type="diverging", robust=True)

    print("Composing animation...")
    save_triptych_gif(
        P_history,
        u_history,
        v_history,
        root / "flow_fields.gif",
        meta=meta,
        p_norm=p_norm,
        u_norm=u_norm,
        v_norm=v_norm,
        times=times,
        fps=12,
        interval=80,
        dpi=170,
        interpolation="nearest",
        transpose=True,
    )
    print(f"{root / 'flow_fields.gif'} written.")
    print("Done.")