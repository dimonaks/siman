import random
import math
import itertools
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # needed for some older Python versions

# -------------------------------------------------------------------------
# 1) Symmetry operations for cubic (BCC) point group: O_h (48 ops)
# -------------------------------------------------------------------------
def generate_o_h_symmetry_ops():
    """
    Generate the 48 symmetry operations (3x3 integer matrices) of the O_h group
    (the full octahedral group). Each matrix has elements in {-1, 0, +1},
    representing all permutations of axes (6) times sign flips (8).
    """
    ops = []
    for perm in itertools.permutations([0,1,2]):
        for signs in itertools.product([+1,-1],[+1,-1],[+1,-1]):
            M = [[0,0,0],[0,0,0],[0,0,0]]
            for row in range(3):
                col = perm[row]
                M[row][col] = signs[row]
            ops.append(M)
    return ops  # 48 operations total

def apply_symmetry(shift, op):
    """
    Apply symmetry operation `op` (3x3) to shift (sx, sy, sz),
    then reduce each coordinate modulo 1 to keep it in [0,1).
    
    shift: (float, float, float)
    op:    3x3 list of lists (integers in {-1, 0, +1})
    returns: tuple (sx', sy', sz') in [0,1).
    """
    s_new = [0.0, 0.0, 0.0]
    for i in range(3):
        s_new[i] = (op[i][0]*shift[0] +
                    op[i][1]*shift[1] +
                    op[i][2]*shift[2])
    # Wrap each coordinate into [0,1)
    for i in range(3):
        s_new[i] %= 1.0
    return tuple(s_new)

def euclid_dist(a, b):
    """Euclidean distance between two 3D points a and b."""
    return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)

def are_sym_equiv_or_close(new_s, existing_s, sym_ops, min_dist):
    """
    Check if 'new_s' is within 'min_dist' of 'existing_s' directly
    OR within 'min_dist' of any of its symmetry-images (via 'sym_ops').
    
    Return True if they are "too close" or "symmetry-equivalent."
    """
    # Direct distance
    if euclid_dist(new_s, existing_s) < min_dist:
        return True
    
    # Check all symmetry images of existing_s
    for op in sym_ops:
        s_image = apply_symmetry(existing_s, op)
        if euclid_dist(new_s, s_image) < min_dist:
            return True
    
    return False


# -------------------------------------------------------------------------
# 2) Generate shifts by subdividing the box, with BCC symmetry + distance checks
# -------------------------------------------------------------------------
def generate_uniform_random_shifts_bcc_sym(m=3, min_dist=0.1, attempt_limit=1000):
    """
    Divide [0,1)^3 into an m x m x m grid of sub-boxes.
    In each sub-box, randomly place exactly one shift in that local region.
    
    Then check:
      - distance to previously accepted shifts (>= min_dist),
      - distance to all their BCC symmetry images (>= min_dist).
    
    If a random candidate is "too close," retry in that sub-box until success
    or until 'attempt_limit' is reached. If no success, we skip that sub-box.
    
    Returns a list of accepted shifts (up to m^3 if all sub-boxes succeed).
    """
    sym_ops = generate_o_h_symmetry_ops()
    shifts = []
    
    # Size of each sub-box along each axis
    sub_size = 1.0 / m
    
    for i in range(m):
        for j in range(m):
            for k in range(m):
                # We'll try up to 'attempt_limit' times to place a shift in this sub-box
                success = False
                for _ in range(attempt_limit):
                    # Generate a random shift in the sub-box (i, j, k)
                    sx = random.uniform(i*sub_size, (i+1)*sub_size)
                    sy = random.uniform(j*sub_size, (j+1)*sub_size)
                    sz = random.uniform(k*sub_size, (k+1)*sub_size)
                    cand = (sx, sy, sz)
                    
                    # Check distance+symmetry vs existing shifts
                    too_close = False
                    for sh in shifts:
                        if are_sym_equiv_or_close(cand, sh, sym_ops, min_dist):
                            too_close = True
                            break
                    
                    if not too_close:
                        # Accept this candidate
                        shifts.append(cand)
                        success = True
                        break
                
                if not success:
                    print(f"Warning: could not place a shift in sub-box {i,j,k} "
                          f"within {attempt_limit} attempts.")
    return shifts


# -------------------------------------------------------------------------
# 3) Plot the resulting shifts in a 3D cube
# -------------------------------------------------------------------------
def plot_shifts_3d(shifts):
    """
    Plot the given list of (sx, sy, sz) shifts in a 3D cube [0,1]^3 using matplotlib.
    """
    xs, ys, zs = zip(*shifts)
    
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    
    ax.scatter(xs, ys, zs)
    
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_zlim(0, 1)
    
    ax.set_xlabel('Shift X')
    ax.set_ylabel('Shift Y')
    ax.set_zlabel('Shift Z')
    
    plt.show()


def equidistant_shifts(mk):
    dk = 1/mk
    shifts = []
    for m1 in range(mk):
        for m2 in range(mk):
            for m3 in range(mk):
                sk1 = dk*m1
                sk2 = dk*m2
                sk3 = dk*m3
                shiftk = [sk1,sk2,sk3]
                shifts.append(shiftk)
                    # print(im, shiftk)

    return shifts


# -------------------------------------------------------------------------
# Example usage
# -------------------------------------------------------------------------
if __name__ == "__main__":
    # Suppose we want a 3x3x3 sub-box division, so up to 27 shifts.
    # We require a min_dist of 0.1 between any two shifts or symmetry images.
    # We allow up to 1000 attempts per sub-box.
    m = 3
    min_dist = 0.11
    attempt_limit = 1000
    
    shifts = generate_uniform_random_shifts_bcc_sym(m=m,
                                                    min_dist=min_dist,
                                                    attempt_limit=attempt_limit)
    print(f"\nGenerated {len(shifts)} shifts from an {m}x{m}x{m} subdivision.\n")
    for i, s in enumerate(shifts, start=1):
        print(f"{i:2d}: {s}")
    
    # Visualize
    plot_shifts_3d(shifts)
