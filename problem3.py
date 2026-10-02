# ---------------------------------------------------------
# TESTING PROBLEM 3
# ---------------------------------------------------------
from helper_functions import local_alignment

answer = local_alignment("pending itch", "unending glitch", lambda x, y: [-1, 1][x == y])
expected = ('ending --itch', 'ending glitch', 9.0)

if (answer == expected):
    print("problem 3 implemented correctly! (hopefully)")
else:
    print(f"expected {expected}, but got {answer}")
