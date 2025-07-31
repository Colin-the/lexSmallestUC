from itertools import product

# sequence = [
#     "2415","3215","3425","1435","2435","1425",
#     "3145","3245","3125","4135","4215","4325",
#     "4315","4235","4125","3415","2345"
# ]

sequence = ["123", "124", "132", "134", "214", "324", "314", "234"]
ROOT = "2345"

def substrings3(x):
    return {
        x[0:3],
        x[1:4],
        x[2:4] + x[0],
        x[3]   + x[0:2]
    }

def substrings2(x):
    return {
        x[0:2],
        x[1:3],
        x[2]   + x[0]
    }

S = { x: substrings2(x) for x in sequence }
print(S)
# 1) build the allowed‐parents lists
allowed = {}
for x in sequence:
    if x == ROOT:
        continue
    allowed[x] = [y for y in sequence if y != x and (S[x] & S[y])]

print(allowed)

invalidCnt = 0
# 2) a quick cycle-&-connectivity checker
def is_valid_tree(parent_map):
    # ensure every node reaches ROOT without revisiting
    for x in sequence:
        seen = set()
        cur = x
        while cur is not None:
            if cur in seen:
                return False   # cycle!
            seen.add(cur)
            cur = parent_map.get(cur, None)
        if ROOT not in seen:
            invalidCnt += 1
            return False       # x didn’t reach root
    return True

# 3) iterate over every choice of parents
solutions = []
nodes = list(allowed.keys())
for picks in product(*(allowed[x] for x in nodes)):
    pm = dict(zip(nodes, picks))
    pm[ROOT] = None

    if not is_valid_tree(pm):
        continue

    # build adjacency
    cand = {x: [] for x in sequence}
    for child, par in pm.items():
        if par is not None:
            cand[par].append(child)
    # ensure the same child‐ordering as in your target list
    for u in cand:
        cand[u].sort(key=lambda c: pos[c])

    # postorder
    out = []
    def post(u):
        for v in cand[u]:
            post(v)
        out.append(u)
    post(ROOT)

    if out == sequence:
        solutions.append(pm)

# 4) report
print("Looked at ",invalidCnt, " invlaid trees")
print(f"Found {len(solutions)} tree(s) matching exactly.")
for sol in solutions:
    print(sol)
