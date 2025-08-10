import matplotlib.pyplot as plt
import networkx as nx
from collections import deque, defaultdict
from typing import Dict, List
import itertools
import math

def plot_tree(
    tree: Dict[str, List[str]],
    root: str,
    filename: str = "tree.png",
    node_size: int = 800,
    font_size: int = 8,
    layer_gap: float = 1.5,
    sibling_gap: float = 1.0,
):
    """
    Draw a hierarchical tree of string nodes, including isolated nodes.

    Parameters
    ----------
    tree : Dict[str, List[str]]
        Mapping from parent node to list of its children.
    root : str
        The node from which to start the hierarchy.  Any nodes not
        reachable from `root` will be laid out at the bottom.
    filename : str
        Where to save the image.
    node_size : int
        Size of each drawn node (in points^2).
    font_size : int
        Font size for node labels.
    layer_gap : float
        Vertical spacing between levels.
    sibling_gap : float
        Minimum horizontal spacing between siblings.
    """
    # --- Build graph and collect all nodes ---
    G = nx.DiGraph()
    all_nodes = set(tree.keys())
    for parent, children in tree.items():
        all_nodes.update(children)
        for child in children:
            G.add_edge(parent, child)
    # include isolated nodes
    for n in all_nodes:
        if n not in G:
            G.add_node(n)

    # --- Compute depth (level) for reachable nodes via BFS from root ---
    depth = {}
    if root in all_nodes:
        depth[root] = 0
        queue = deque([root])
        while queue:
            u = queue.popleft()
            for v in tree.get(u, []):
                if v not in depth:
                    depth[v] = depth[u] + 1
                    queue.append(v)

    # Any node not reached from root, give it depth = max_depth + 1
    max_depth = max(depth.values()) if depth else 0
    for n in all_nodes:
        if n not in depth:
            depth[n] = max_depth + 1

    # --- Group nodes by level ---
    levels = defaultdict(list)
    for n, d in depth.items():
        levels[d].append(n)
    n_levels = len(levels)

    # --- Assign positions ---
    pos = {}
    # we'll center each level on x=0, spread siblings evenly
    for d in range(n_levels):
        group = sorted(levels[d])
        count = len(group)
        # total width for this level
        width = max((count - 1) * sibling_gap, 0)
        start = -width / 2
        for i, node in enumerate(group):
            x = start + i * sibling_gap
            y = -d * layer_gap
            pos[node] = (x, y)

    # --- Draw ---
    plt.figure(figsize=(max(8, n_levels * 1.5), max(6, len(all_nodes) * 0.2)))
    nx.draw(
        G,
        pos,
        with_labels=True,
        labels={n: n for n in all_nodes},
        node_size=node_size,
        font_size=font_size,
        arrowsize=12,
        arrowstyle='<|-',
    )
    plt.title("Tree")
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    plt.show()



# nodes = [
#     "2315",
#     "2415","3215","3425","1435","2435","1425",
#     "3145","3245","3125","4135","4215","4325",
#     "4315","4235","4125","3415","2345"
# ]

nodes = [
    "1234", "1235", "1243", "1245", "1324", "1325", "1342", "1345", "2135", "2143", "2145", "2314", "2315",
    "2415","3215","3425","1435","2435","1425",
    "3145","3245","3125","4135","4215","4325",
    "4315","4235","4125","3415","2345"
]

ROOT = "2345"

def inOrder(node):
    if node is None:
        return False
    pos = 1
    for num in node:
        if int(num) > pos + 1:
            #print(num,"is out of order in",node)
            return False
        pos+= 1
    return True

def missing(node, n):
    for num in range(1, n + 1):
        if str(num) not in node:
            return num
    return -1

def parent(node: str, n):
    # The root has no parent or it's parent is it's self
    if node == ROOT:
        return ROOT
    if node == "2415":
        return "2315"
    if node == "3145":
        return "3125"
    
    # Determine the missing number
    missingNum = missing(node, n)
    # If there was no missing number then we have a error 
    # in the permutation and can't find it's parent
    if missingNum == -1:
        return None
    
    if inOrder(node):
        candidates = []
        currentVal = int(node)
        # Try inserting the missing digit at every possible index
        for i in range(len(node)):
            cand = node[:i] + str(missingNum) + node[i+1:len(node)]
            #print(cand)
            # Only consider those greater than the current node numerically
            if int(cand) > currentVal:
                candidates.append(cand)

        # Return the lexicographically smallest new string that is greater than
        # the current permutation
        if candidates:
            return min(candidates, key=lambda x: int(x))
        return None
    else:
        checkInv = False
        # IF the current perm does not contain 1 and is out of order
        if missingNum == n:
            return None
        elif missingNum == n - 1:  
            toReplace = 1
        elif missingNum == math.ceil(n/2):
            toReplace = n + 1
            checkInv = True
        else:
            toReplace = missingNum + 1

        # Try replacing m+1 with m where m is the missing symbol
        for i in range(len(node)):
            if int(node[i]) == toReplace:
                par = node[:i] + str(missingNum) + node[i+1:len(node)]

        if missingNum == 1 or missingNum == 2:
            return par
        
        # if checkInv:
        #     # Generate the other possible parent by replacing 1 with m
        #     # Try replacing m+1 with m where m is the missing symbol
        #     for i in range(len(node)):
        #         if int(node[i]) == 1:
        #             par2 = node[:i] + str(missingNum) + node[i+1:len(node)]

        #     # Now we want par or par2 depending on...?
        
        short = par[:len(par)-1]
        print(short)
        valid = True
        for i in range(1,len(short) - 1):
            if int(short[i]) < int(short[i+1]):
                continue
            else:
                valid = False

        if valid:
            return par
        

        for i in range(len(node)):
            if int(node[i]) == 1:
                par = node[:i] + str(missingNum) + node[i+1:len(node)]

        # If we are transtioning between out and in order
        # then we need to make sure that this is the right
        # time to do so
        if inOrder(par):
            # Strip out the last symbol
            short = par[:len(par)-1]

            # Now we want the lex greatest IN order combonation
            # of these numbers 
            validInversion = True

            # numb = [int(s) for s in list(short)]
            # numb.sort()
            # combs = itertools.permutations(numb, len(short))

            # #print(numb)
            # for combo in combs:
            #     #print(combo," and int ",int("".join(map(str, combo))))
            #     if int("".join(map(str, combo))) > int(short) and inOrder(combo):
            #         #print("\t\tNot vaild!")
            #         validInversion = False
            #         break
            """
            for i in range(len(short)-1):
                # compare adj number
                if short[i] > short[i+1]:
                    continue
                # If the next number in the sequence is greater than the current one
                else:
                    # Build the new perm that will be lex greater then the current one
                    new = short[:max(i,0)]
                    new += short[i+1]
                    new += short[i]
                    new += short[i+2:]

                    print(new)
                    # If this new perm is in order then we can't swtich from out of
                    # order to in order and we need to appily the next rule to find the parent
                    if inOrder(new):
                        print(new,"is not in order")
                        validInversion = False
                    else:
                        continue
            """
            

            #print(short)
            if validInversion:
                return par
            
            return None
        return par


if __name__ == '__main__':
    test_cases = [
        ("124", 4),
        ("2315", 5),
        ("1234", 5),  
        ("2415", 5)   # inOrder check fails, so None
    ]
    # for node, n in test_cases:
    #     print(f"parent({node}, {n}) = {parent(node, n)}")

    for node in nodes:
        p = parent(node, 5)
        #if p is not None and not inOrder(node):
        if not inOrder(node):
            if inOrder(p):
                print("\tswtich! ",node,"has parent of", p)
            else:
                if p in nodes:
                    print(node,"has parent of",p)
                else:
                    print("Parent out of",node,"out of set!",p)


    tree: Dict[str, List[str]] = defaultdict(list)

    # For every node, compute its parent and add it as a child
    for node in nodes:
        p = parent(node, 5)
        # skip the root mapping to self, or include if you want a self‐loop
        if p and p != node:
            tree[p].append(node)

    # Make sure every node appears as a key, even if it has no children
    for node in nodes + [ROOT]:
        tree.setdefault(node, [])

    # for parent_node, kids in tree.items():
    #     print(f"{parent_node!r} → {kids}")

    plot_tree(tree, root="2345", filename="my_tree.png")



# 4325, 4315, 4215, 4235, 4135, 4125, 3125, 3425, 3415, 2415