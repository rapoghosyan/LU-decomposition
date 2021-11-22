import time


print('Вводите матрицу по строкам:')
a = list(map(float, input().split()))
n = len(a)
A = []
LR = []
for i in range(n):
    A.append([])
    LR.append([])
    for j in range(n):
        A[i].append(0)
        LR[i].append(0)
for j in range(n):
    A[0][j] = a[j]
for i in range(1, n):
    a = list(map(float, input().split()))
    for j in range(n):
        A[i][j] = a[j]
r = pow(0.1, 7)
# LR алгоритм нахождения собственных значений для произвольной матрицы
start_time = time.time()
sd = 0
while True:
    # Для матрицы Ak, найдем его LR разложение, и норми R - Ak, L - E
    #                                      (||A|| = (∑a[i][j]^2)^1/2 )
    ok = 1
    error1 = 0
    error2 = 0
    LR[0][0] = A[0][0]
    for i in range(1, n):
        LR[0][i] = A[0][i]
        if A[i][0] == 0:
            LR[i][0] = 0
        elif LR[0][0] == 0:
            ok = 0
            break
        else:
            LR[i][0] = A[i][0] / LR[0][0]
            error2 += LR[i][0] * LR[i][0]
    if ok == 1:
        for i in range(1, n):
            for l in range(i, n):
                s = 0
                for j in range(0, i):
                    s += (LR[i][j] * LR[j][l])
                LR[i][l] = A[i][l] - s
                error1 += (s * s)
            for l in range(i + 1, n):
                s = 0
                for j in range(0, i):
                    s += (LR[l][j] * LR[j][i])
                if A[l][i] - s == 0:
                    LR[l][i] = 0
                elif LR[i][i] == 0:
                    ok = 0
                    break
                else:
                    LR[l][i] = (A[l][i] - s) / LR[i][i]
                    error2 += LR[l][i] * LR[l][i]
            if ok == 0:
                break
    if ok == 1:
        for i in range(1, n):
            for j in range(i):
                error1 += (A[i][j] * A[i][j])
        error1 = pow(error1, 0.5)
        error2 = pow(error2, 0.5)
        if error1 < r and error2 < r:
            print('Собственные значения:')
            for i in range(n):
                print(' λ', i + 1, ' = ', A[i][i], sep='')
            break
        else:
            # Переходим на следующую матрицу - Ak+1 = Rk * Lk
            for i in range(n):
                for j in range(n):
                    s = 0
                    for l in range(max(i, j), n):
                        if j == l:
                            s += LR[i][l]
                        else:
                            s += LR[i][l] * LR[l][j]
                    A[i][j] = s
            # Если было сделано сдвыг, то Ak+1 = Rk * Lk + λE, где λ → 0
            if sd == 1:
                sd = 0
                for i in range(n):
                    A[i][i] += r
    else:
        # Переходим на следующую матрицу - Ak+1 = Ak - λE, где λ → 0
        sd = 1
        for i in range(n):
            A[i][i] -= r
end_time = time.time()
print('\n', 'Время работы - ', "%.10f" % (end_time - start_time), ' с.',  sep='')
