package williamfiset;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.HashMap;

public class combinations {
    // This method generates all bit sets of size n where r bits
    // are set to one. The result is returned as a list of integer masks.
    public static List<Integer> combinationsF(int r, int n) {
        List<Integer> subsets = new ArrayList<>();
        combinationsHelper(0, 0, r, n, subsets);
        return subsets;
    }

    // To find all the combinations of size r we need to recurse until we have
    // selected r elements (aka r = 0), otherwise if r != 0 then we still need to
    // select
    // an element which is found after the position of our last selected element
    private static void combinationsHelper(int set, int at, int r, int n, List<Integer> subsets) {

        // Return early if there are more elements left to select than what is
        // available.
        int elementsLeftToPick = n - at;
        if (elementsLeftToPick < r)
            return;

        // We selected 'r' elements so we found a valid subset!
        if (r == 0) {
            subsets.add(set);
        } else {
            for (int i = at; i < n; i++) {
                // Try including this element
                set ^= (1 << i);

                combinationsHelper(set, i + 1, r - 1, n, subsets);

                // Backtrack and try the instance where we did not include this element
                set ^= (1 << i);
            }
        }
    }

    public static void main(String[] args) {

        List<Integer> arrI = combinationsF(3, 10);

        for (Integer element : arrI) {
            System.out.println(Integer.toBinaryString(element));
        }

    }
}
