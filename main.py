def input_equations():
    print("=== 3x3 연립방정식 Ax = b 입력 ===")

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
