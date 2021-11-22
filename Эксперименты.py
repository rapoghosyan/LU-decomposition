import time
import random
import sys


def sob_zn(a, sd):
    global mmn, r
    n = len(a)
    if n == 2:
        d = pow(a[0][0] + a[1][1], 2) - 4 * (a[0][0] * a[1][1] - a[0][1] * a[1][0])
        if d < 0:
            x = (a[0][0] + a[1][1]) / 2
            return [str(x) + " + " + str(pow(-d, 0.5) / 2) + " i", str(x) + " - " + str(pow(-d, 0.5) / 2) + " i"]
        elif d == 0:
            return [(a[0][0] + a[1][1]) / 2, (a[0][0] + a[1][1]) / 2]
        else:
            x1 = (a[0][0] + a[1][1] - pow(d, 0.5)) / 2
            x2 = (a[0][0] + a[1][1] + pow(d, 0.5)) / 2
            if abs(x1) < abs(x2):
                x1, x2 = x2, x1
            return [x1, x2]
    if n == 1:
        return a[0]
    if sd == 0:
        for i in range(n - 2, -1, -1):
            if abs(a[i + 1][i]) <= r * mmn:
                temp1 = []
                for j in range(i + 1):
                    temp1.append([])
                    for l in range(i + 1):
                        temp1[j].append(a[j][l])
                temp2 = []
                for j in range(n - i - 1):
                    temp2.append([])
                    for l in range(i + 1, n):
                        temp2[j].append(a[j + i + 1][l])
                return sob_zn(temp1, 0) + sob_zn(temp2, 0)
        i = n - 1
        while i >= 0 and a[i][i] == 0:
            i -= 1
        if i == -1:
            sd = r
        else:
            sd = a[i][i]
        for i in range(n):
            a[i][i] -= sd
    lr = []
    for i in range(n):
        lr.append([])
        for j in range(n):
            lr[i].append(0)
    ok = 1
    lr[0][0] = a[0][0]
    for i in range(1, n):
        lr[0][i] = a[0][i]
    for i in range(1, n):
        if a[i][i - 1] == 0:
            lr[i][i - 1] = 0
        elif lr[i - 1][i - 1] == 0:
            ok = 0
            break
        else:
            lr[i][i - 1] = a[i][i - 1] / lr[i - 1][i - 1]
        for l in range(i, n):
            s = lr[i][i - 1] * lr[i - 1][l]
            lr[i][l] = a[i][l] - s
    if ok == 1:
        for l in range(n - 1):
            a[0][l] = lr[0][l] + lr[0][l + 1] * lr[l + 1][l]
        a[0][n - 1] = lr[0][n - 1]
        for i in range(1, n):
            for l in range(i, n - 1):
                a[i][l] = lr[i][l] + lr[i][l + 1] * lr[l + 1][l]
            a[i][n - 1] = lr[i][n - 1]
            a[i][i - 1] = lr[i][i] * lr[i][i - 1]
        for i in range(n):
            a[i][i] += sd
        return sob_zn(a, 0)
    else:
        for i in range(n):
            a[i][i] -= r
        return sob_zn(a, sd + r)


sys.setrecursionlimit(1000000)
kal_pav = 100
for test in range(100, 101, 1):
    print(test)
    test_time = 0
    min_test = 10000000
    max_test = 0
    pochti_tre = 0
    for pav in range(kal_pav):
        nx = test
        matrix = []
        for ix in range(nx):
            matrix.append([])
            for jx in range(nx):
                matrix[ix].append(0)
        for ix in range(nx):
            for jx in range(ix, nx):
                matrix[ix][jx] = random.uniform(-1000, 1000)
        r = pow(0.1, 7)
        start_time = time.time()
        ax = []
        for jx in range(1, nx):
            ax.append(0)
        for ix in range(nx - 2):
            for jx in range(1, nx - ix - 1):
                for kx in range(nx - ix - 1):
                    ax[kx] = matrix[kx + ix + 1][ix]
                x0 = pow(ax[0] * ax[0] + ax[jx] * ax[jx], 0.5)
                if x0 != 0:
                    cox = ax[0] / x0
                    six = -ax[jx] / x0
                    # matrix = Tij * matrix
                    for kx in range(nx):
                        xi = matrix[ix + 1][kx]
                        xj = matrix[jx + ix + 1][kx]
                        matrix[ix + 1][kx] = xi * cox - xj * six
                        matrix[jx + ix + 1][kx] = xi * six + xj * cox
                    # matrix = matrix * Tji
                    for kx in range(nx):
                        xi = matrix[kx][ix + 1]
                        xj = matrix[kx][jx + ix + 1]
                        matrix[kx][ix + 1] = xi * cox - xj * six
                        matrix[kx][jx + ix + 1] = xi * six + xj * cox
        pochti_tre += time.time() - start_time
        mmn = 0
        for ix in range(nx):
            sx = 0
            for jx in range(max(0, ix - 1), nx):
                sx += abs(matrix[ix][jx])
            mmn = max(mmn, sx)
        sob_zn(matrix, 0)
        end_time = time.time()
        min_test = min(min_test, end_time - start_time)
        max_test = max(max_test, end_time - start_time)
        test_time += end_time - start_time
    print("%.17f" % (min_test))
    print("%.17f" % (max_test))
    print("%.17f" % (test_time / kal_pav))
    print("%.17f" % (pochti_tre / kal_pav))
