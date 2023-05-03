m = [[i+j*10 for j in range(10)] for i in range(10)]
i = 3
j = 8
m[i][j] = 'X'
for row in m:
    print(row)

k = [x for x in range(i+1, j-1)]
print(k)
I = [m[i][x] for x in k]
J = [m[x+1][j] for x in k]
print(I)
print(J)
print([x + y for x, y in zip(I, J)])



