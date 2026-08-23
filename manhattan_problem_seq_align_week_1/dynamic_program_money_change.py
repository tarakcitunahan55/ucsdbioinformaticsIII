def dp_change(money, coins):
    """
    Solve the Change Problem using dynamic programming (bottom-up).

    Unlike recursive_change, which starts at 'money' and works DOWN
    to 0 (recomputing overlapping subproblems many times), this version
    starts at 0 and works UP to 'money', filling in a table as it goes.
    Each subproblem is solved exactly once and reused instantly afterward.

    Parameters:
        money - the amount of change needed
        coins - list of available coin denominations in decreasing order

    Returns:
        Minimum number of coins needed to make 'money'.
    """

    # min_num_coins[m] will store the minimum coins needed to make amount m.
    # We need entries for every amount from 0 up to 'money', so the list has (money + 1) slots. Start all of them as infinity ("not yet known").
    min_num_coins = [float('inf')] * (money + 1)

    # Base case: it costs 0 coins to make 0 change. This is the ONE value we know for certain before any computation -- everything else builds on it.
    min_num_coins[0] = 0

    # Build the table up from m=1 to m=money. By the time we're computing min_num_coins[m], every smaller amount (m - coin) has ALREADY been fully computed and stored 
    # (no recursion, no re-deriving old answers, just a direct lookup.)
    for m in range(1, money + 1):
        # For this amount m, try every coin and see which one gives the cheapest way to make m.
        for coin in coins:
            if m >= coin:
                # min_num_coins[m - coin] is already known (computed earlier in this same loop, since m - coin < m). +1 accounts for using this coin itself.
                if min_num_coins[m - coin] + 1 < min_num_coins[m]:
                    min_num_coins[m] = min_num_coins[m - coin] + 1

    # After the loop finishes, min_num_coins[money] holds our final answer.
    return min_num_coins[money]


if __name__ == "__main__":
    coins = [20, 18, 5, 3, 1]
    amount = 17583


    result = dp_change(amount, coins)
    print(f"Minimum coins to make {amount} using {coins}: {result}")