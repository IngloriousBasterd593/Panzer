n = int(input("Enter the value of n: "))

lines = n
stars = 1

for i in range(n):
    # Print spaces
    for j in range(lines):
        print(" ", end="")

    # Print stars
    for k in range(stars):
        print("*", end="")

    print()  # Move to the next line

    stars += 2
    lines -= 1
