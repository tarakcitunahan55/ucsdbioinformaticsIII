def dp_change_space_efficient(money, coins):
    """
    Same logic as DPChange, but the array size is capped at max(coins)
    instead of growing to size (money + 1).

    Why this works: computing MinNumCoins[m] only ever needs values from
    m - coin, for coin in coins. Since every coin <= max(coins), the oldest
    value we could possibly need is at index (m - max(coins)). Anything
    older than that is never referenced again -- so we can safely overwrite
    it. We use (index % max_coin) to "wrap around" and reuse old slots.
    """
    max_coin = max(coins)

    # Only need max_coin slots -- not money+1 slots.
    min_num_coins = [float('inf')] * max_coin
    min_num_coins[0] = 0   # base case: 0 coins for amount 0, stored at slot 0

    for m in range(1, money + 1):
        slot = m % max_coin

        # IMPORTANT: reset this slot before computing it. It currently holds
        # a stale value left over from (m - max_coin) iterations ago --
        # we're about to overwrite it with the answer for the NEW amount m.
        min_num_coins[slot] = float('inf')

        for coin in coins:
            if m >= coin:
                prev_slot = (m - coin) % max_coin
                if min_num_coins[prev_slot] + 1 < min_num_coins[slot]:
                    min_num_coins[slot] = min_num_coins[prev_slot] + 1

    return min_num_coins[money % max_coin]