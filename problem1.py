# ---------------------------------------------------------
# TESTING PROBLEM 1
# ---------------------------------------------------------
from helper_functions import global_alignment

answer = global_alignment("abracadabra", "dabarakadara", lambda x, y: [-1, 1][x == y])
expected = ('-ab-racadabra', 'dabarakada-ra', 5.0)

if (answer == expected):
    print("problem 1 implemented correctly! (hopefully)")
else:
    print(f"expected {expected}, but got {answer}")
