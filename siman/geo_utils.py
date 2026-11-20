import numpy as np

from siman.geo import image_distance, xred2xcart
from siman.header import printlog


def _endpoint_cart(st):
    """Cartesian positions of a structure, in angstroms.

    Uses st.xcart when the structure already has Cartesian coordinates.
    Otherwise converts the reduced coordinates st.xred with the cell st.rprimd.
    """
    return st.xcart if st.xcart else xred2xcart(st.xred, st.rprimd)


def _ordering_failure(st1, st2, still_tol):
    """None when the atom order is valid. Otherwise the reason it is not."""
    if st1.natom != st2.natom:
        return 'the endpoints have %d and %d atoms' % (st1.natom, st2.natom)

    elements1 = st1.get_elements()
    elements2 = st2.get_elements()
    mismatched = [
        (i, elements1[i], elements2[i])
        for i in range(st1.natom) if elements1[i] != elements2[i]
    ]
    if mismatched:
        shown = ', '.join('%d: %s/%s' % item for item in mismatched[:8])
        return 'element order differs at ' + shown

    for vector1, vector2 in zip(st1.rprimd, st2.rprimd):
        if np.linalg.norm(np.array(vector1) - np.array(vector2)) > 0.05:
            return 'the endpoint cells differ'

    cart1 = _endpoint_cart(st1)
    cart2 = _endpoint_cart(st2)
    distances = [
        image_distance(x1, x2, st1.rprimd)[0]
        for x1, x2 in zip(cart1, cart2)
    ]
    moving = [i for i, distance in enumerate(distances) if distance > still_tol]
    if len(moving) != 1:
        largest = sorted(range(st1.natom), key=lambda i: distances[i], reverse=True)[:5]
        detail = ', '.join(
            '%d %s %.3f A' % (i, elements1[i], distances[i]) for i in largest)
        return (
            'expected one atom beyond %.2f A, found %d. Largest shifts: %s'
            % (still_tol, len(moving), detail))
    return None


def _ordering_summary(st1, st2, still_tol):
    """One-line description of a valid endpoint order.

    Names the atom that moves farthest between st1 and st2, its displacement
    in angstroms, and the still_tol limit that the other atoms stay inside.
    """
    elements = st1.get_elements()
    cart1 = _endpoint_cart(st1)
    cart2 = _endpoint_cart(st2)
    distances = [
        image_distance(x1, x2, st1.rprimd)[0]
        for x1, x2 in zip(cart1, cart2)
    ]
    atom = max(range(st1.natom), key=lambda i: distances[i])
    return (
        'atom %d (%s) moves %.3f A; the other atoms stay within %.2f A'
        % (atom, elements[atom], distances[atom], still_tol))


def _assign_columns(cost):
    """Match each row of a square cost matrix to one unused column.

    The returned list has one entry per row: assignment[row] is the chosen
    column. Each column is used once. The pairing is the one with the lowest
    total cost.

    SciPy's linear_sum_assignment solves this exactly when SciPy is installed.
    Otherwise the cheapest unused row-column pairs are taken in order, which
    can miss the optimal pairing. _reorder_endpoint uses this to match atoms
    of one element by distance.
    """
    try:
        from scipy.optimize import linear_sum_assignment
    except ImportError:
        linear_sum_assignment = None
    rows = len(cost)
    if linear_sum_assignment is not None:
        chosen_rows, chosen_cols = linear_sum_assignment(np.array(cost, float))
        assignment = [None] * rows
        for row, col in zip(chosen_rows, chosen_cols):
            assignment[row] = col
        return assignment
    pairs = sorted(
        (cost[i][j], i, j) for i in range(rows) for j in range(rows))
    used_rows, used_cols = set(), set()
    assignment = [None] * rows
    for _, row, col in pairs:
        if row in used_rows or col in used_cols:
            continue
        assignment[row] = col
        used_rows.add(row)
        used_cols.add(col)
    return assignment


def _reorder_endpoint(st1, st2):
    """Permute st2 so each atom matches the same element in st1 by position."""
    elements1 = st1.get_elements()
    elements2 = st2.get_elements()
    if sorted(elements1) != sorted(elements2):
        raise RuntimeError('reorder cannot fix different element counts')
    cart1 = _endpoint_cart(st1)
    cart2 = _endpoint_cart(st2)
    new_order = [-1] * st1.natom
    for element in sorted(set(elements1)):
        index1 = [i for i, name in enumerate(elements1) if name == element]
        index2 = [i for i, name in enumerate(elements2) if name == element]
        cost = [
            [image_distance(cart1[i], cart2[j], st1.rprimd)[0] for j in index2]
            for i in index1
        ]
        assignment = _assign_columns(cost)
        for row, col in enumerate(assignment):
            new_order[index1[row]] = index2[col]
    if sorted(new_order) != list(range(st2.natom)):
        raise RuntimeError('reorder did not assign every atom')
    return st2.reorder(new_order)


def _reorder_failure(st1, st2, host_tol=5.0):
    """None when one atom migrates and every other atom is within host_tol."""
    elements1 = st1.get_elements()
    if elements1 != st2.get_elements():
        mismatched = [
            (i, elements1[i], st2.get_elements()[i])
            for i in range(st1.natom) if elements1[i] != st2.get_elements()[i]
        ]
        shown = ', '.join('%d: %s/%s' % item for item in mismatched[:8])
        return 'element order differs at ' + shown
    cart1 = _endpoint_cart(st1)
    cart2 = _endpoint_cart(st2)
    distances = [
        image_distance(x1, x2, st1.rprimd)[0]
        for x1, x2 in zip(cart1, cart2)
    ]
    migrant = max(range(st1.natom), key=lambda i: distances[i])
    tied = [
        i for i, distance in enumerate(distances)
        if abs(distance - distances[migrant]) < 1e-6
    ]
    if len(tied) != 1:
        return 'more than one atom shares the largest shift'
    hosts = [distances[i] for i in range(st1.natom) if i != migrant]
    if any(distance > host_tol for distance in hosts):
        largest = sorted(range(st1.natom), key=lambda i: distances[i], reverse=True)[:5]
        detail = ', '.join(
            '%d %s %.3f A' % (i, elements1[i], distances[i]) for i in largest)
        return 'atoms other than the migrant exceed %.2f A: %s' % (host_tol, detail)
    return None


def order_neb_endpoints(st1, st2, still_tol=1.0, reorder_tol=5.0):
    """Check the given order. On failure, reorder st2 and allow 5 A for the host.

    Both results are written to the log. The returned structure is st2, or
    st2 permuted onto st1 when the second check passes.
    """
    original = _ordering_failure(st1, st2, still_tol)
    if original is None:
        passed = _ordering_summary(st1, st2, still_tol)
        printlog('check_ordering: original order passed: ' + passed, imp='y')
        return st2

    printlog('check_ordering: original order failed: ' + original, imp='y')
    try:
        reordered = _reorder_endpoint(st1, st2)
        reordered_problem = _reorder_failure(st1, reordered, reorder_tol)
    except RuntimeError as error:
        reordered_problem = str(error)
        reordered = None
    if reordered_problem is None:
        passed = _ordering_summary(st1, reordered, reorder_tol)
        printlog('check_ordering: reordered order passed: ' + passed, imp='y')
        return reordered

    printlog('check_ordering: reordered order failed: ' + reordered_problem, imp='y')
    raise RuntimeError(
        'check_ordering failed. Original: %s. After reorder: %s'
        % (original, reordered_problem))


def check_neb_ordering(st1, st2, still_tol=1.0):
    """Require the same atom order, with one migrating atom.

    Corresponding atoms must be the same element. Every atom except one must
    stay within *still_tol* angstroms under the minimum-image convention.
    """
    problem = _ordering_failure(st1, st2, still_tol)
    if problem is not None:
        raise RuntimeError('check_ordering: ' + problem)
    printlog('check_ordering: ' + _ordering_summary(st1, st2, still_tol), imp='y')
    return st2
