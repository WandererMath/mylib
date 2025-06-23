import os
import bisect


def get_paths_ends_with_something(folder_path, suffix):
    return [os.path.abspath(os.path.join(folder_path, f)) for f in os.listdir(folder_path) if f.endswith(suffix)]


def closest_binary_search(nums_sorted, x):
    pos = bisect.bisect_left(nums_sorted, x)
    if pos == 0:
        return nums_sorted[0]
    if pos == len(nums_sorted):
        return nums_sorted[-1]
    before = nums_sorted[pos - 1]
    after = nums_sorted[pos]
    if after - x < x - before:
        return after
    else:
        return before

if __name__ == "__main__":
    # Example usage:
    nums = [10, 22, 14, 3, 76, 54]
    nums = sorted(nums)
    

    queries = [5, 50, 1, 100]
    for q in queries:
        print(closest_binary_search(nums, q))
        break

def _auto_grouping(samples, key= lambda x: x):
    result={}
    for s in samples:
        if key(s)[0] not in result:
            result[key(s)[0]]=[s]
        else:
            result[key(s)[0]].append(s)
    return result