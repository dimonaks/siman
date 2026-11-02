import re
from collections import defaultdict

def parse_bands(outcar="OUTCAR", tol=1e-6):
    data, spin, k = {}, 1, None
    rk = re.compile(r"k-point\s+(\d+)")
    rb = re.compile(r"^\s*(\d+)\s+([-+.\dEe]+)\s+([-+.\dEe]+)")

    with open(outcar, errors="ignore") as f:
        for line in f:
            if "spin component" in line:
                spin = int(line.split()[-1])
            elif m := rk.search(line):
                k = int(m.group(1))
            elif k is not None and (m := rb.match(line)):
                b, e, occ = int(m[1]), float(m[2]), float(m[3])
                data[spin, k, b] = e, occ

    by_spin = defaultdict(lambda: defaultdict(list))
    for (s, kp, b), (e, occ) in data.items():
        by_spin[s][b].append((e, occ))

    for s, bands in sorted(by_spin.items()):
        occupied, empty, mixed = [], [], []
        occ_e, empty_e = [], []

        for b, values in sorted(bands.items()):
            states = [o > tol for _, o in values]
            energies = [e for e, _ in values]

            if all(states):
                occupied.append(b)
            elif not any(states):
                empty.append(b)
            else:
                mixed.append(b)

            occ_e += [e for e, o in values if o > tol]
            empty_e += [e for e, o in values if o <= tol]

        print(f"\nSpin {s}: {len(bands)} bands, "
              f"{len({kp for ss, kp, _ in data if ss == s})} k-points")
        print(f"  Occupied bands : {len(occupied):4d}")
        print(f"  Empty bands    : {len(empty):4d}")
        print(f"  Mixed bands    : {len(mixed):4d}")
        print(f"  Occupied range : {min(occ_e):7.1f} ... {max(occ_e):7.1f} eV")
        print(f"  Empty range    : {min(empty_e):7.1f} ... {max(empty_e):7.1f} eV")