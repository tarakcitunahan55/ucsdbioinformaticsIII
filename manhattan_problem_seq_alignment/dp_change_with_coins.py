def dp_change_with_coins(money, coins):
    """
    Computes the minimum number of coins AND reconstructs the actual
    coins used, by keeping a second array that remembers which coin
    was responsible for each amount's best answer.
    """
    min_num_coins = [float('inf')] * (money + 1)
    min_num_coins[0] = 0

    # last_coin_used[m] = which coin denomination gave the best answer for m.
    # This is the extra bookkeeping needed to reconstruct the actual coins.
    last_coin_used = [None] * (money + 1)

    for m in range(1, money + 1):
        for coin in coins:
            if m >= coin:
                if min_num_coins[m - coin] + 1 < min_num_coins[m]:
                    min_num_coins[m] = min_num_coins[m - coin] + 1
                    last_coin_used[m] = coin   # remember the winning coin for this m

    # Backtrack from 'money' down to 0, following the trail of coins used.
    coins_used = []
    remaining = money
    while remaining > 0:
        coin = last_coin_used[remaining]
        coins_used.append(coin)
        remaining -= coin   # step back to the amount before this coin was added

    return min_num_coins[money], coins_used


# ---- Example usage ----
if __name__ == "__main__":
    coins = [25, 10, 5, 1]
    amount = 37

    count, used = dp_change_with_coins(amount, coins)
    print(f"Minimum coins for {amount}: {count}")
    print(f"Coins used: {used}  (sum = {sum(used)})")