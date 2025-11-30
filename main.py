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
    
def input_equations():
    print("3x3 연립방정식 Ax = b 입력")

    A = []
    for i in range(3):
        row = list(map(float, input(f"A의 {i+1}번째 행 (3개 입력): ").split()))
        A.append(row)

    b = list(map(float, input("b 벡터 (3개 입력): ").split()))

    return A, b


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
