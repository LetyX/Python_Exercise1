import numpy as np


def get_evals(A):
    e, evecs = np.linalg.eig(A)
    e.sort()
    return e


def count_repeats(v):
    e = []
    degen = []
    x = round(v[0], 3)
    k = 1
    for i in range(1, len(v)):
        if abs(v[i] - x) <= 0.001:
            k += 1
        else:
            e.append(x)
            degen.append(k)
            x = round(v[i], 3)
            k = 1
    e.append(x)
    degen.append(k)
    print("The energy levels and their degeneracies are:\n")
    for j in range(0, len(e)):
        print(str(e[j].real).rjust(7), degen[j], sep='   ')


def lpoly(n):
    hm = np.zeros((n, n))
    hm[0, 1] = -1
    hm[n - 1, n - 2] = -1
    for i in range(1, n - 1):
        hm[i, i - 1] = -1
        hm[i, i + 1] = -1
    return hm


def cpoly(n):
    hm = np.zeros((n, n))
    for i in range(0, n):
        hm[i, (i - 1) % n] = -1
        hm[i, (i + 1) % n] = -1
    return hm


def plats(f):
    # set up connectivity list
    if f == 4:  # tetrahedron
        v = 4
        conn = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]]
    elif f == 6:  # cube
        v = 8
        conn = [[1, 3, 7], [0, 2, 4], [1, 3, 5], [0, 2, 6], [1, 5, 7], [2, 4, 6], [3, 5, 7], [0, 4, 6]]
    else:  # dodecahedron
        v = 20
        conn = [[1, 4, 5], [0, 2, 6], [1, 3, 7], [2, 4, 8], [0, 3, 9], [0, 10, 14], [1, 10, 11], [2, 11, 12],
                [3, 12, 13], [4, 13, 14], [5, 6, 15], [6, 7, 16], [7, 8, 17], [8, 9, 18], [5, 9, 19],
                [10, 16, 19], [11, 15, 17], [12, 16, 18], [13, 17, 19], [14, 15, 18]]

    # Build Huckel matrix for the Platonic solid
    hm = np.zeros((v, v))
    for i in range(0, v):
        for j in range(i + 1, v):
            if j in conn[i]:
                hm[i, j] = -1
                hm[j, i] = -1

    return hm


def bucky():
    # Build lists showing connectivity for 2 halves of the buckminsterfullerene
    conn1 = []
    conn2 = []
    for i in range(0, 5):  # Build a pentagonal graph
        conn1.append([(i - 1) % 5, (i + 1) % 5])
        conn2.append([59 - (i - 1) % 5, 59 - (i + 1) % 5])

    for i in range(5, 10):  # Attach one new vertex to each of the pentagon vertices
        conn1[i % 5].append(i)
        conn1.append([i % 5])
        conn2[i % 5].append(59 - i)
        conn2.append([59 - i % 5])

    for i in range(10, 20):
        # Add a new pair of vertices between each of the previously added vertices to form 5 hexagonal faces
        if i % 2 == 0:
            conn1[int(i / 2)].append(i)
            conn1.append([int(i / 2), i + 1])
            conn2[int(i / 2)].append(59 - i)
            conn2.append([59 - int(i / 2), 59 - (i + 1)])
        else:
            conn1[int((i + 1) / 2) % 5 + 5].append(i)
            conn1.append([i - 1, int((i + 1) / 2) % 5 + 5])
            conn2[int((i + 1) / 2) % 5 + 5].append(59 - i)
            conn2.append([59 - (i - 1), 59 - (int((i + 1) / 2) % 5 + 5)])

    for i in range(20, 30):  # Add 5 new pairs of vertices to form 5 pentagonal faces
        conn1[((i + 1) % 10) + 10].append(i)
        conn2[((i + 1) % 10) + 10].append(59 - i)
        if i % 2 == 0:
            conn1.append([((i + 1) % 10) + 10, i + 1])
            conn2.append([59 - (((i + 1) % 10) + 10), 59 - (i + 1)])
        else:
            conn1.append([((i + 1) % 10) + 10, i - 1])
            conn2.append([59 - (((i + 1) % 10) + 10), 59 - (i - 1)])

    # Join the two halves
    for i in range(0, 10):
        conn1[20 + i % 10].append(30 + (i - 1) % 10)
        conn2[20 + i % 10].append(20 + (i - 1) % 10)

    conn2.reverse()

    # Build Huckel matrix
    hm = np.zeros((60, 60))
    for i in range(0, 30):
        for j in range(i + 1, 60):
            if j in conn1[i]:
                hm[i, j] = -1
                hm[j, i] = -1

    for i in range(0, 30):
        for j in range(i + 31, 60):
            if j in conn2[i]:
                hm[30 + i, j] = -1
                hm[j, 30 + i] = -1

    return hm


print("Choose a type of molecule:")
print("1 linear polyene", "2 cyclic polyene", "3 sp2 hybridised Platonic solid", "4 Buckminsterfullerene", sep=', ')
option = int(input("Enter a number: "))
while option not in range(1, 5):
    print("Your option must be an integer between 1 and 4 inclusive.")
    option = int(input("Enter a number: "))

if option == 1 or option == 2:
    y = int(input("Enter the number of atoms: "))
    while y <= 0:
        print("The number of atoms must be a non-zero positive integer.")
        y = int(input("Enter the number of atoms: "))
    if option == 1:
        h = lpoly(y)
    else:
        h = cpoly(y)
elif option == 3:
    y = int(input("Enter the number of faces: "))
    while y not in [4, 6, 12]:
        print("The number of faces has to be 4 (tetrahedron), 6 (cube) or 12 (dodecahedron).")
        y = int(input("Enter the number of faces: "))
    h = plats(y)
else:
    h = bucky()

count_repeats(get_evals(h))
