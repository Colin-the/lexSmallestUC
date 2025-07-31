# ——————————————————————————————
# 1) define your sequence & root
# ——————————————————————————————
sequence = [
    "2415","3215","3425","1435","2435","1425",
    "3145","3245","3125","4135","4215","4325",
    "4315","4235","4125","3415","2345"
]
ROOT = "2345"

# ——————————————————————————————
# 2) compute the 4 rotated length‐3 substrings
# ——————————————————————————————
def substrings3(x):
    return {
        x[0:3],
        x[1:4],
        x[2:4] + x[0],
        x[3]   + x[0:2]
    }

S = { x: substrings3(x) for x in sequence }

print(S)
# ——————————————————————————————
# 3) the parent rule as described
# ——————————————————————————————
def parent(x):
    if x == ROOT:
        return None
    best = None
    best_key = (-1, None)  # (common‐substr‐as‐int, -int(y)) to maximize
    for y in sequence:
        if y == x:
            continue
        common = S[x] & S[y]
        if not common:
            continue
        m = max(int(s) for s in common)
        key = (m, -int(y))
        if key > best_key:
            best_key = key
            best = y
    return best


def possibleParents(x):
    if x == ROOT:
        return None
    related = []
    for y in sequence:
        if y == x:
            continue
        common = S[x] & S[y]

        # If they don't share 3 symbols they can't be related
        if not common:
            continue

        related.append(y)
    return related

# ——————————————————————————————
# 4) build the adjacency list
# ——————————————————————————————
tree = { x: [] for x in sequence }
for x in sequence:
    p = parent(x)
    poss = possibleParents(x)
    print("Parent of ",x," is ",p)
    print(x," could be related to: ",poss)
    if p is not None:
        tree[p].append(x)

print(tree)
# now sort each node's children in the order they appear in `sequence`,
# so that our postorder will match your exact list:
pos = { x:i for i,x in enumerate(sequence) }
for x in tree:
    tree[x].sort(key=lambda c: pos[c])

# ——————————————————————————————
# 5) postorder traversal & check
# ——————————————————————————————
def postorder(u, out):
    for v in tree[u]:
        postorder(v, out)
    out.append(u)

out = []
postorder(ROOT, out)

print(out)
assert out == sequence, "❌ postorder did *not* match!"
print("✅ postorder matches exactly!")



# {'2415': ['2435'], '3215': ['3245'], '3425': ['3415'], '1435': ['1425'], '2435': ['2415'], '1425': ['1435'], '3145': ['3125'], '3245': ['3215'], '3125': ['3145'], '4135': ['4125'], '4215': ['4235'], '4325': ['4315'], '4315': ['4325'], '4235': ['4215'], '4125': ['4135'], '3415': ['3425'], '2345': []}