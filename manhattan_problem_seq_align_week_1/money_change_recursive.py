def recursive_change(money, coins):
    """
    Recursively compute the minimum number of coins needed to make
    change for a given amount of money, using the given coin denominations.

    Parameters:
        money - the amount of change needed
        coins - a list of available coin denominations, e.g. [25,10,5,1] in decreasing order

    Returns:
        The minimum number of coins needed to make up 'money'.

    WARNING: this direct recursive translation is extremely slow for
    larger amounts, because it recomputes the same subproblems over
    and over (exponential time)- a dynamic programming version fixes this.
    """

    # Base case: making change for 0 requires 0 coins.
    if money == 0:
        return 0

    # Start with "infinity" as a placeholder -- any real solution found in the loop will be smaller than this.
    min_num_coins = float('inf')

    # Try using each coin denomination as the "last" coin used.
    for i in range(len(coins)):
        # Only consider this coin if it doesn't exceed the amount left.
        if money >= coins[i]:
            # Recursively solve for the remaining amount after using this coin.
            num_coins = recursive_change(money - coins[i], coins)

            # If using this coin leads to a better (smaller) solution, keep it. (+1 accounts for the coin we just used.)
            if num_coins + 1 < min_num_coins:
                min_num_coins = num_coins + 1

    return min_num_coins


if __name__ == "__main__":
    coins = [25,10,5,1] 
    amount = 37

    result = recursive_change(amount, coins)
    print(f"Minimum coins to make {amount} using {coins}: {result}")