def input_equations():
    print("3x3 연립방정식 Ax = b 입력")

    A = []
    for i in range(3):
        row = list(map(float, input(f"A의 {i+1}번째 행 (3개 입력): ").split()))
        A.append(row)

    b = list(map(float, input("b 벡터 (3개 입력): ").split()))

    return A, b

def gaussianElimination(A, b):
  
  
    M = [row[:] + [bx] for row, bx in zip(A, b)]
    rows, cols = 3, 4
    pivot_row = 0

    # 1. 전진 소거
    for j in range(3): # 열 순회 (변수 x1, x2, x3)
        if pivot_row >= rows: break

        # 현재 열에서 가장 큰 절댓값 원소를 피벗으로 선택
        max_val, max_row = 0, -1
        for i in range(pivot_row, rows):
            if abs(M[i][j]) > max_val:
                max_val = abs(M[i][j])
                max_row = i
        
        # 피벗이 0에 가까우면 다음 열로 이동
        if max_val < 1e-9: continue 

        # 행 교환
        M[pivot_row], M[max_row] = M[max_row], M[pivot_row]

        # 소거
        pivot = M[pivot_row][j]
        for i in range(pivot_row + 1, rows):
            factor = M[i][j] / pivot
            for k in range(j, cols):
                M[i][k] -= factor * M[pivot_row][k]
                if abs(M[i][k]) < 1e-9: M[i][k] = 0.0 # 부동 소수점 오차 처리
        
        pivot_row += 1

    # 2. 후진 소거
    # 피벗을 1로 만들고 피벗 위의 원소들을 0으로 소거
    
    pivots = []
    for i in range(rows):
        for j in range(3):
            if abs(M[i][j]) > 1e-9:
                pivots.append((i, j))
                
                # 피벗을 1로 만들기
                pivot_val = M[i][j]
                for k in range(j, cols):
                    M[i][k] /= pivot_val
                
                break
    
    # 후진하며 피벗 위의 원소 소거
    for i, j in reversed(pivots):
        for row_idx in range(i):
            factor = M[row_idx][j]
            for k in range(j, cols):
                M[row_idx][k] -= factor * M[i][k]
                if abs(M[row_idx][k]) < 1e-9: M[row_idx][k] = 0.0

    return M

P = [[0.0]*3 for _ in range(3)]
L = [[0.0]*3 for _ in range(3)]
U = [[0.0]*3 for _ in range(3)]

P = [[0.0]*3 for _ in range(3)]
L = [[0.0]*3 for _ in range(3)]
U = [[0.0]*3 for _ in range(3)]

def lu_decomposition(A):
    global P, L, U

    # P를 단위행렬로 초기화, L과 U는 0으로 초기화; U는 A의 복사본으로 시작
    P = [[0.0]*3 for _ in range(3)]
    for i in range(3):
        P[i][i] = 1.0
    L = [[0.0]*3 for _ in range(3)]
    U = [list(map(float, row)) for row in A]  # A를 float 형태로 복사

    # LU 분해, 부분 피벗팅
    for k in range(3):
        # k번째 열에서 k행부터 2행까지 절댓값 최대 피벗 행 찾기
        pivot_row = max(range(k, 3), key=lambda r: abs(U[r][k]))

        # 피벗 행이 현재 k행이 아니면, U와 P의 행을 교환하고
        # 이미 계산된 L의 이전 열 값도 교환
        if pivot_row != k:
            U[k], U[pivot_row] = U[pivot_row], U[k]
            P[k], P[pivot_row] = P[pivot_row], P[k]
            # 이미 계산된 곱셈 계수(L) 교환 (열 0..k-1)
            for j in range(k):
                L[k][j], L[pivot_row][j] = L[pivot_row][j], L[k][j]

        # L의 대각 원소는 1로 설정
        L[k][k] = 1.0

        # 아래 행 제거를 위한 곱셈 계수 계산 및 U 행 갱신
        for i in range(k+1, 3):
            L[i][k] = U[i][k] / U[k][k]
            for j in range(k, 3):
                U[i][j] -= L[i][k] * U[k][j]

    # 글로벌 변수에 다시 저장 (지역 변수로 이미 계산했지만 참조 보장)
    globals()['P'] = P
    globals()['L'] = L
    globals()['U'] = U
    # L,U 출력
    print("L =")
    for row in L:
        print("  [" + ", ".join(f"{x: .6f}" for x in row) + "]")
    print()

    print("U =")
    for row in U:
        print("  [" + ", ".join(f"{x: .6f}" for x in row) + "]")
    print()

def back_substitution(U, y):
    x = [0.0] * 3
    for i in range(2, -1, -1):   # 2,1,0 순서
        sum_value = y[i]
        for j in range(i+1, 3):  # 이미 구한 x
            sum_value -= U[i][j] * x[j]
        x[i] = sum_value / U[i][i]
    return x

def forward_substitution(L, b):
    y = [0.0] * 3
    for i in range(3):
        sum_value = b[i]
        for j in range(i):     # 이미 구한 y
            sum_value -= L[i][j] * y[j]
        y[i] = sum_value / L[i][i]
    return y

def main():
    A, b = input_equations()

    sol_type, info = check_solution_type(A, b)

    if sol_type == "NO_SOLUTION":
        print("해가 없습니다.")
    
    elif sol_type == "UNIQUE":
        L, U = lu_decomposition(A)
        y = forward_substitution(L, b)
        x = back_substitution(U, y)
        print("유일해:", x)

    else:  # INFINITE
        print("무수히 많은 해:", info)
