nodes = [
    "1234", "1235", "1243", "1245", "1324", "1325", "1342", "1345", "2135", "2143", "2145", "2314", "2315",
    "2415","3215","3425","1435","2435","1425",
    "3145","3245","3125","4135","4215","4325",
    "4315","4235","4125","3415","2345"
]
ROOT = "2345"

def inOrder(node):
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
    if node == ROOT:
        return None
    
    # Determine the missing number
    missingNum = missing(node, n)
    if missingNum == -1:
        return None
    
    if inOrder(node):
        # 12341234 (for 1234)
        candidates = []
        current_val = int(node)
        # Try inserting the missing digit at every possible index
        for i in range(len(node)):
            cand = node[:i] + str(missingNum) + node[i+1:len(node)]
            #print(cand)
            # Only consider those greater than the current node numerically
            if int(cand) > current_val:
                candidates.append(cand)

        #print(candidates)
        # Return the lexicographically smallest by numeric value
        if candidates:
            return min(candidates, key=lambda x: int(x))
        return None
    else:
        return None

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
        if parent(node, 5) is not None:
            print(node,"has parent")
