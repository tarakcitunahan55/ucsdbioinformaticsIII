def hanoi_towers(n, start_peg, destination_peg, move_count=None):
    """
    Recursively solve the Tower of Hanoi for n disks.

    Parameters:
        n               - number of disks to move
        start_peg       - the peg the disks start on (1, 2, or 3)
        destination_peg - the peg the disks need to end up on (1, 2, or 3)
        move_count      - a mutable list used to count total moves
    
    Why move_count[0] += 1 only in base case and step 2?
    Step 1 and Step 3 are not disk moves themselves — they're recursive calls that eventually contain disk moves (via their own base cases and their own Step 2's, several layers down). If we incremented the counter at the point of the recursive call itself, we'd be double-counting: once when the call is made, and again for every move that call generates internally.
    every recursive call eventually bottoms out at a base case, and each base case contributes exactly one increment
   
    For example, tracing n=2:

hanoi_towers(2, 1, 3) calls hanoi_towers(1, 1, 2) → this hits the base case, increments once (that's the "Step 1" work for n=2, but the increment happens inside the n=1 call, not in the n=2 call's Step 1 line).
Back in hanoi_towers(2, 1, 3), Step 2 runs: increments once.
Then hanoi_towers(1, 2, 3) → base case again, increments once.

Total: 3 increments for 3 actual moves — correct for 2 disks (2² - 1 = 3).

So the pattern is: every recursive call eventually bottoms out at a base case, and each base case contributes exactly one increment. Step 2 in the middle of each call contributes one more.
     """

    # move_count lets us track total moves across all recursive calls.
    # We use a list (mutable) instead of an int (immutable) -> immutable problem in recursion (normally would work fine in non-recursion) since each recursive call gets its own private copy of move_count (starting fresh at whatever value was passed in). When that nested call increments it, that increment only exists inside that call's own local scope
    if move_count is None:
        move_count = [0]

    # ---- BASE CASE ----
    # If there's just 1 disk, there's nothing to "clear out of the way" — just move it directly from start to destination.
    if n == 1:
        move_count[0] += 1
        print(f"Move {move_count[0]}: Move disk 1 from peg {start_peg} to peg {destination_peg}")
        return move_count[0]

    # ---- RECURSIVE CASE ----
    # Figure out which peg is the "helper" (transit) peg.
    # It's whichever peg ISN'T start or destination.
    transit_peg = 6 - start_peg - destination_peg

    # Step 1: Move the top (n-1) disks out of the way onto the transit peg.
    # From this sub-call's point of view, destination_peg is just a helper — it isn't disks n-1's final peg yet.
    hanoi_towers(n - 1, start_peg, transit_peg, move_count)

    # Step 2: Move the single largest disk (disk n) from start to destination.
    # This is the disk we were "clearing space" for in Step 1.
    move_count[0] += 1
    print(f"Move {move_count[0]}: Move disk {n} from peg {start_peg} to peg {destination_peg}")

    # Step 3: Move the (n-1) disks from the transit peg onto the destination peg,
    # placing them on top of the disk we just moved.
    # Here, start_peg becomes the helper peg for this sub-call.
    hanoi_towers(n - 1, transit_peg, destination_peg, move_count)

    return move_count[0]


# ---- Example usage ----
if __name__ == "__main__":
    num_disks = 4
    print(f"Solving Tower of Hanoi for {num_disks} disks (peg 1 -> peg 3):\n")

    total_moves = hanoi_towers(num_disks, start_peg=1, destination_peg=3)

    print(f"\nTotal moves used: {total_moves}")
    print(f"Theoretical minimum (2^n - 1): {2**num_disks - 1}")