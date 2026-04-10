import argparse
import numpy as np
from scipy.io import loadmat


def main():
    p = argparse.ArgumentParser(
        description="Export one array from a MATLAB .mat file to a compressed .npz (key arr_0)."
    )
    p.add_argument("mat_path", help="Input .mat file")
    p.add_argument("npz_path", help="Output .npz file")
    p.add_argument(
        "--var",
        default="arr_0",
        help="Variable name inside the .mat file (default: arr_0)",
    )
    args = p.parse_args()

    d = loadmat(args.mat_path)
    if args.var not in d:
        keys = [k for k in d if not k.startswith("__")]
        raise SystemExit(
            f"variable {args.var!r} not in {args.mat_path}; available: {keys}"
        )
    arr = np.asarray(d[args.var])
    np.savez_compressed(args.npz_path, arr_0=arr)


if __name__ == "__main__":
    main()
