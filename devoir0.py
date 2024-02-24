def solve_ratio(phi):
    # résoudre le système suivant :
    # x1 + x2 = 1
    # phi * x1 + (1-phi) * x2 = 1
    if phi == 1 / 2:
        return None
    x2 = (1 - phi) / (1 - 2 * phi)
    x1 = 1 - x2
    print("Verif :", x1 + x2, (phi * x1) + (1 - phi) * x2)
    return x1, x2


print(solve_ratio(1))
print(solve_ratio(2))
print(solve_ratio(3))
print(solve_ratio(4))
print(solve_ratio(5))
print(solve_ratio((1 + 5**0.5) / 2))
