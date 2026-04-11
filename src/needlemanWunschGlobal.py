def needleman_wunsch(seq1, seq2):
    # Настройки очков
    match = 1
    mismatch = -1
    gap = -2

    # Определяем размеры будущей таблицы
    n = len(seq1)
    m = len(seq2)


    # Создаем матрицу, заполненную нулями
    # Initialize the scoring matrix with zeros
    score_matrix = [[0 for _ in range(m + 1)] for _ in range(n + 1)]

    # Fill the first column and row with gap penalties
    # Заполняем первый столбец
    for i in range(n + 1):
        score_matrix[i][0] = i * gap

    # Заполняем первую строку
    for j in range(m + 1):
        score_matrix[0][j] = j * gap


    # Fill the matrix using dynamic programming
    # Заполняем остальную часть матрицы
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # Считаем очки для трёх вариантов:
            # Calculate scores for match/mismatch and gaps
            if seq1[i-1] == seq2[j-1]:
                score = match
            else:
                score = mismatch

            # Диагональ (совпадение/несовпадение)
            diagonal = score_matrix[i-1][j-1] + score
            # Сверху (гэп в seq1)
            up = score_matrix[i-1][j] + gap
            # Слева (гэп в seq2)
            left = score_matrix[i][j-1] + gap

            # Выбираем максимум из трёх
            # Store the highest score in the current cell
            score_matrix[i][j] = max(diagonal, up, left)

            # Обратный путь (Traceback)
            # Traceback to find the optimal alignment
    align1 = ""
    align2 = ""
    i = n
    j = m

    while i > 0 or j > 0:
        current_score = score_matrix[i][j]

        # Check if the move came from the diagonal
        # 1. Проверяем, пришли ли мы по диагонали (совпадение/несовпадение)
        if i > 0 and j > 0 and (current_score == score_matrix[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)):
            align1 += seq1[i-1]
            align2 += seq2[j-1]
            i -= 1
            j -= 1

            # Check if the move came from above (gap in seq2)
        # 2. Проверяем, пришли ли мы сверху (гэп в seq2)
        elif i > 0 and current_score == score_matrix[i-1][j] + gap:
            align1 += seq1[i-1]
            align2 += "-"
            i -= 1
            # Otherwise, the move came from the left (gap in seq1)
        # 3. Иначе мы пришли слева (гэп в seq1)
        else:
            align1 += "-"
            align2 += seq2[j-1]
            j -= 1

    # Переворачиваем строки, так как мы шли с конца
    # Reverse the strings since we traced back from the end
    print(f"\nAlignment Result:")
    print(align1[::-1])
    print(align2[::-1])

    print(f"Final Score: {score_matrix[n][m]}")

    return align1[::-1], align2[::-1], score_matrix[n][m]

if __name__ == "__main__":
    # Test the algorithm with example sequences
    result_matrix = needleman_wunsch("GATTACA", "GATCA")

    print("Final Score Matrix:")
    for row in result_matrix:
      print(row)

