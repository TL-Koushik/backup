import java.util.*;

class sorts {
  public static long Pow(int X, int N) {
    if (N <= 1) {
      return X;
    }
    return (long) X * (long) Math.pow(X, N - 1);
  }

  public static String removeConsecutiveDuplicates(String s) {
    String se = "" + s.charAt(0);
    if (se.charAt(se.length() - 1) != s.charAt(0)) {
      se = se.substring(se.length() - 1) + s.charAt(0);
    }
    System.out.println(se);
    if (s.length() == 1) {
      return se;
    }
    removeConsecutiveDuplicates(s.substring(1, s.length()));
    return se;
  }

  public static int[] selection(int[] arr) {
    int temp = 0;
    for (int i = 0; i < arr.length; i++) {
      for (int j = i; j < arr.length; j++) {
        if (arr[i] > arr[j]) {
          temp = arr[j];
          arr[j] = arr[i];
          arr[i] = temp;
        }
      }
    }
    return arr;
  }

  public static int reverse(int n) {
    if (n % 10 == n) {
      return n;
    }
    return n % 10 * (int) Math.pow(10, (int) (Math.log10(n))) + reverse(n / 10);
  }

  public static boolean helper(int i, int[] a) {
    if (i == a.length - 1) {
      return true;
    }
    return a[0] <= a[1] && helper(++i, a);
  }

  public static boolean checksorted(int[] a) {
    if (a.length == 1) {
      return true;
    }
    if (a[a.length - 2] > a[a.length - 1]) {
      return false;
    }
    return checksorted(Arrays.copyOf(a, a.length - 1));
  }

  public static int linearsearch(int[] a, int key, int idx) {
    if (idx == a.length - 1) {
      return -1;
    }
    if (key == a[idx]) {
      return idx;
    } else {
      return linearsearch(a, key, ++idx);
    }
  }

  public static ArrayList linearsearchbyalloccurance(int[] a, int key, int idx) {
    ArrayList<Integer> ls = new ArrayList<>();
    if (idx == a.length) {
      return ls;
    }
    if (a[idx] == key) {
      ls.add(idx);
    }
    return linearsearchbyalloccurance(a, key, idx - 1);
  }

  public static int rotatedsortedarray(int[] a, int key) {
    int s = 0, l = a.length - 1;
    int mid = 0;
    int pivot = a[a.length - 1];
    while (s <= l) {
      mid = s + (l - s) / 2;
      if (a[mid] > a[s]) {
        s = mid + 1;
        pivot = mid;
      }

    }
    if (a[a.length - 1] < key) {
      s = 0;
      l = pivot;
    } else {
      s = pivot + 1;
      l = a.length - 1;
    }
    while (s <= l) {
      mid = (s + l) / 2;
      if (a[mid] == key) {
        return mid;
      } else if (a[mid] < key) {
        s = mid + 1;
      } else {
        l = mid - 1;
      }

    }
    return -1;
  }

  public static int rotatedsortedarraybyrec(int[] a, int key, int s, int e) {
    if (s > e) {
      return -1;
    }
    int mid = (s + (e - s)) / 2;
    if (a[mid] == key) {
      return mid;
    }
    if (key >= a[s]) {
      if (key <= a[mid]) {
        return rotatedsortedarraybyrec(a, key, s, mid - 1);
      } else {
        return rotatedsortedarraybyrec(a, key, mid + 1, e);
      }
    }
    if (key >= a[mid] && key <= a[e]) {
      return rotatedsortedarraybyrec(a, key, mid + 1, e);// arr=.key=2,3
    }
    return rotatedsortedarraybyrec(a, key, s, mid - 1);
  }

  public static int[] bubblesortbyrecursion(int[] a, int d, int s) {
    if (s == 0) {
      return a;
    }
    if (s == d) {
      return bubblesortbyrecursion(a, 0, s - 1);
    }
    if (a[d] > a[d + 1]) {
      int temp = a[d];
      a[d] = a[d + 1];
      a[d + 1] = temp;
    }
    return bubblesortbyrecursion(a, ++d, s);

  }

  public static int[] selectionbyrecursion(int[] a, int search, int indx, int temp) {
    if (search == a.length - 1) {
      return a;
    }
    if (indx == a.length) {
      indx = a[temp];
      a[temp] = a[search];
      a[search] = indx;
      search++;
      return selectionbyrecursion(a, search, search + 1, temp);// 3,4,5,6,7,8,9,1,2,3
    }
    if (a[temp] > a[indx]) {
      temp = indx;
    }
    return selectionbyrecursion(a, search, ++indx, temp);

  }

  public static int[] twoSum(int[] a, int target) {
    HashMap<Integer, Integer> s = new HashMap<>();
    int[] r = new int[2];
    for (int i = 0; i < a.length; i++) {
      if (s.containsKey(a[i])) {
        r[0] = 1 + (s.get(target - a[i]));
        r[1] = ++i;
        return r;
      }
      s.put(a[i], i);
      s.put(target - a[i], -1);
    }
    return r;
  }

  static ArrayList<ArrayList<Integer>> l = new ArrayList<>();

  public static void threeSum(int[] nums) {
    Arrays.sort(nums);
    int j, k;
    for (int i = 0; i < nums.length && nums[i] <= 0; i++) {
      k = nums.length - 1;
      j = i + k / 2;
      if (i == k || j == k || i == j) {
        return;
      }
      while (k < nums.length) {
        if (j == i) {
          k--;
          j++;
        }
        if (j == k) {
          break;
        }
        if (nums[i] + nums[j] + nums[k] == 0) {
          System.out.println(+nums[i] + " " + nums[j] + " " + nums[k]);
          break;
        } else if (nums[i] + nums[j] + nums[k] < 0) {
          j++;
        } else {
          j--;
        }
      }
    }

  }

  public static boolean isPalindrome(String s) {
    if (s.length() <= 1) {
      return true;
    }
    s = s.toLowerCase();
    char l = s.charAt(0), r = s.charAt(s.length() - 1);
    if (Character.isLetterOrDigit(l) && Character.isLetterOrDigit(r)) {
      if (l != r) {
        return false;
      }
      s = s.substring(1, s.length() - 1);
    } else if (!Character.isLetterOrDigit(l)) {
      s = s.substring(1, s.length());
    } else if (!Character.isLetterOrDigit(r)) {
      s = s.substring(0, s.length() - 1);
    }
    return isPalindrome(s);
  }

  public static int[] combine(int[] left, int[] right) {
    int[] arr = new int[left.length + right.length];
    int k = 0, l = 0, r = 0;
    while (k < arr.length && r < right.length && l < left.length) {
      if (right[r] < left[l]) {
        arr[k++] = right[r++];
      } else {
        arr[k++] = left[l++];
      }
    }
    while (r < right.length) {
      arr[k++] = right[r++];
    }
    while (l < left.length) {
      arr[k++] = left[l++];
    }
    return arr;

  }

  public static int[] mergesort(int[] a) {
    if (a.length <= 1) {
      return a;
    }
    int m = a.length / 2;
    int[] left = mergesort(Arrays.copyOfRange(a, 0, m));
    int[] rigth = mergesort(Arrays.copyOfRange(a, m, a.length));
    return combine(left, rigth);
  }

  public static String backspaceCompare(String s, String t) {
    String sres = "" + s.charAt(0);
    int i = 1;
    while (i < s.length()) {
      if (s.charAt(i) == '#') {
        sres = sres.substring(0, sres.length() - 1);
      } else {
        sres = sres + s.charAt(i);
      }
      i++;
    }
    return sres;

  }

  public static int[] getSecondOrderElements(int n, int[] a) {
    int[] re = new int[2];
    int max = 0;
    int secondleast = Integer.MAX_VALUE;
    int min = Integer.MAX_VALUE;
    for (int i = 0; i < n; i++) {
      if (a[i] < min) {
        re[1] = min;
        min = a[i];
      } else if (a[i] < re[1] && a[i] != min) {
        re[1] = a[i];
      }
      if (max < a[i]) {
        re[0] = max;
        max = a[i];
      }
    }
    return re;

  }

  public static int getlongestsubarray(int[] a, int k) {
    int sum = Integer.MIN_VALUE;
    int max = 0;
    for (int i = 0; i < a.length; i++) {
      for (int j = i; j < a.length; j++) {
        sum = sum + a[j];
        if (sum == k) {
          max = Math.max(max, j - i + 1);
        }
      }
    }
    return max;
  }

  public static int getLongestSubarray(int[] nums, int k) {
    // optimal if the array consists negative and positives
    int sum = 0;
    int max = 0;
    int remining = 0;
    HashMap<Integer, Integer> previoussum = new HashMap<>();
    for (int i = 0; i < nums.length; i++) {
      sum += nums[i];
      if (sum == k) {
        max = Math.max(max, i + 1);
      }
      remining = k - sum;
      if (previoussum.containsKey(remining)) {
        max = Math.max(max, i - previoussum.get(remining));
      }
      if (!previoussum.containsKey(sum)) {
        previoussum.put(sum, i);
      }
    }
    return max;
  }

  public static int getLongestSubarraypostive(int[] nums, int k) {
    int i = 0, j = 0, sum = 0, max = 0;
    while (j < nums.length) {
      sum += nums[j];
      while (sum > k && i <= j) {
        sum -= nums[i];
        i++;
      }
      if (sum == k) {
        max = Math.max(max, j - i + 1);
      }
      j++;
    }
    return max;
  }

  public static String read(int[] book, int target) {
    HashMap<Integer, Integer> remining = new HashMap<>();
    for (int i = 0; i < book.length; i++) {
      if (remining.containsKey(book[i])) {
        return "YES";
      }

      remining.put(target - book[i], i);

    }
    return "NO";
  }

  public static void sort012(int[] a) {
    int i = 0;
    int j = 0;
    int k = a.length - 1;
    int temp = 0;
    while (j <= k) {
      if (a[j] == 0) {
        temp = a[i];
        a[i] = a[j];
        a[j] = temp;
        i++;
        j++;
      } else if (a[j] == 1) {
        j++;
      }

      else {
        temp = a[k];
        a[k] = a[j];
        a[j] = temp;
        k--;
      }

    }
    for (i = 0; i < a.length; i++) {
      System.out.println(a[i]);
    }
  }

  public static int majorityElement(int[] v) {
    // moores algorithm
    int count = 0;
    int element = 0;
    for (int i = 0; i < v.length; i++) {
      if (count == 0) {
        element = v[i];
        count++;
      } else if (v[i] == element) {
        count++;
      } else {
        count--;
      }
    }
    count = 0;
    for (int i = 0; i < v.length; i++) {
      if (element == v[i]) {
        count++;
      }
    }
    if ((v.length / 2) < count) {
      return element;
    }
    return -1;

  }

  public static int maxSubArray(int[] nums) {// kadane algorithm
    int max = Integer.MIN_VALUE;
    int sum = 0;
    int s = 0, end = 0, ansstart = 0;
    for (int i = 0; i < nums.length; i++) {
      if (sum == 0) {
        s = i;
      }
      sum += nums[i];
      if (max < sum) {
        max = sum;
        ansstart = s;
        end = i;
      }
      if (sum < 0) {
        sum = 0;
      }
    }
    for (int i = ansstart; i <= end; i++) {
      System.out.println(nums[i]);
    }
    return max;

  }

  public static int[] rearrangeArray(int[] nums) {
    int i = 0;
    int j = 0;
    int temp = 0;

    while (i < nums.length) {
      if (i % 2 == 1 && nums[i] > 0) {
        j = i + 1;
        while (nums[j] > 0 && j < nums.length) {
          j++;
        }
        temp = nums[i];
        nums[i] = nums[j];
        nums[j] = temp;
      } else if (i % 2 == 0 && nums[i] < 0) {
        j = i + 1;
        while (nums[j] < 0 && j < nums.length) {
          j++;
        }
        temp = nums[i];
        nums[i] = nums[j];
        nums[j] = temp;
      } else {
        i++;
      }
    }
    return nums;

  }

  public static int[] reverse(int[] nums, int start) {
    int end = nums.length - 1;
    while (start <= end) {
      int temp = nums[start];
      nums[start] = nums[end];
      nums[end] = temp;
      start++;
      end--;
    }
    return nums;
  }

  public static int[] nextPermutation(int[] nums) {
    int i = nums.length - 2;
    int lower = 0;
    while (i >= 0) {
      if (nums[i] < nums[i + 1]) {
        lower = i;
        break;

      }
      i--;
    }
    if (i == -1) {
      nums = reverse(nums, 0);
      return nums;
    }
    i = nums.length - 1;
    while (nums[lower] < nums[i]) {
      i--;
    }
    int temp = nums[lower];
    nums[lower] = nums[i];
    nums[i] = temp;
    return reverse(nums, lower);

  }

  public static int longestSuccessiveElements(int[] a) {
    Arrays.sort(a);
    int count = 1;
    int seq = 1;
    for (int i = 0; i < a.length - 1; i++) {
      if (a[i] == a[i + 1]) {
        continue;
      }
      if (a[i] + 1 == a[i + 1]) {
        count++;
      } else {
        count = 1;
      }
      seq = Math.max(seq, count);
    }
    return seq;
  }

  // public static List< Integer > superiorElements(int []a) {
  // ArrayList<Integer> ele=new ArrayList<>();
  // ele.add(a[a.length-1]);
  // int lastgreater=a[a.length-1];
  // for(int i=a.length-2;i>=0;i--){
  // if(a[i]>lastgreater){
  // ele.add(a[i]);
  // lastgreater=a[i];
  // }

  // }
  // return ele;

  // }

  public static int[][] zeroMatrix(int[][] matrix, Integer n, Integer m) {
    int col = 1;
    for (int i = 0; i < matrix[0].length; i++) {
      for (int j = 0; j < matrix.length; j++) {
        if (matrix[i][j] == 0) {
          matrix[i][0] = 0;
          if (j != 0) {
            matrix[0][j] = 0;
          } else {
            col = 0;
          }
        }
      }
    }
    for (int i = 1; i < matrix[0].length; i++) {
      for (int j = 1; j < matrix.length; j++) {
        if (matrix[i][j] != 0) {
          if (matrix[i][0] == 0 || matrix[0][j] == 0) {
            matrix[i][j] = 0;
          }
        }
      }
    }
    if (matrix[0][0] == 0) {
      for (int i = 0; i < matrix[0].length; i++) {
        matrix[0][i] = 0;
      }
    }
    if (col == 0) {
      for (int i = 0; i < matrix.length; i++) {
        matrix[i][0] = 0;
      }
    }
    return matrix;

  }

  public static int[][] anti90matrix(int[][] m) {
    int[][] c = new int[m.length][m[0].length];
    for (int i = 2; i >= 0; i--) {
      for (int j = 0; j < 3; j++) {
        c[i][j] = m[j][2 - i];
      }
    }
    return c;
  }

  public static int[][] degreematspaceop(int[][] matrix) {
    int temp = 0;
    for (int i = 0; i < matrix[0].length; i++) {
      for (int j = i + 1; j < matrix.length; j++) {
        temp = matrix[i][j];
        matrix[i][j] = matrix[j][i];
        matrix[j][i] = temp;

      }
    }
    int s = 0;
    int e = matrix.length - 1;
    for (int i = 0; i < matrix[0].length; i++) {
      while (s <= e) {
        temp = matrix[i][s];
        matrix[i][s++] = matrix[i][e];
        matrix[i][e--] = temp;
      }
      s = 0;
      e = matrix.length - 1;
    }
    return matrix;
  }

  public static List<List<Integer>> generate(int numRows) {
    List<List<Integer>> l = new ArrayList<>();
    for (int i = 0; i < numRows; i++) {
      ArrayList<Integer> row = new ArrayList<>();
      for (int j = 0; j <= i; j++) {
        if (j == 0 || j == i) {
          row.add(1);
        } else {
          int num = l.get(i - 1).get(j - 1) + l.get(i - 1).get(j);
          row.add(num);
        }
      }
      l.add(row);
    }
    return l;

  }

  public static int[] spiralmatrix(int[][] m) {
    int left = 0, top = 0, bottom = m.length - 1, rigth = m[0].length - 1;
    int index = 0;
    int[] list = new int[m.length * m[0].length];
    while (left <= rigth && top <= bottom) {
      for (int i = left; i <= rigth; i++) {
        System.out.println(m[top][i]);
        list[index] = m[top][i];
        index++;
      }
      top++;
      for (int i = top; i <= bottom; i++) {
        System.out.println(m[i][rigth]);
        list[index] = m[i][rigth];
        index++;
      }
      rigth--;
      if (top <= bottom) {
        for (int i = rigth; i >= left; i--) {
          System.out.println(m[bottom][i]);
          list[index] = m[bottom][i];
          index++;
        }
        bottom--;
      }
      if (left <= rigth) {
        for (int i = bottom; i >= top; i--) {
          System.out.println(m[i][left]);
          list[index] = m[i][left];
          index++;
        }
        left++;
      }
    }
    return list;

  }

  public static void triplet(int n, int[] arr) {
    int j = 0;
    int k = 0;
    Arrays.sort(arr);
    ArrayList<ArrayList<Integer>> tri = new ArrayList<>();
    for (int i = 0; i < arr.length; i++) {
      if (i != 0 && arr[i] == arr[i - 1])
        continue;
      j = i + 1;
      k = arr.length - 1;
      while (j < k) {
        if ((arr[i] + arr[j] + arr[k]) < 0) {
          j++;
        } else if ((arr[i] + arr[j] + arr[k]) > 0) {
          k--;
        } else {
          // ArrayList<Integer> now=Arrays.asList(arr[i],arr[j],arr[k]);
          // tri.add(now);
          System.out.print(" " + arr[i] + " " + arr[j] + " " + arr[k]);
          j++;
          k--;
          while (arr[j] == arr[j - 1] && j < k)
            j++;
          while (arr[k] == arr[k + 1] && j < k)
            k--;
        }
      }
    }
  }

  public static int subarraysWithSumK(int[] a, int k) {
    // Write your code here
    int count = 0;
    int xor = 0;
    HashMap<Integer, Integer> hs = new HashMap<>();
    hs.put(0, 1);
    for (int i = 0; i < a.length; i++) {
      xor ^= a[i];// this xor all the elements like adding all the elements into the xor var
      int x = xor ^ k;// this is used to get the element which is stoping to get the subarray of k by
                      // this the xor can filter the elements which is stoping from making the array
                      // the k
      if (hs.containsKey(x)) {
        count += hs.get(x);
      }
      if (hs.containsKey(xor)) {
        hs.put(xor, hs.get(xor) + 1);
      } else {
        hs.put(xor, 1);
      }
    }
    return count;
  }

  // public static List< List< Integer > > mergeOverlappingIntervals(int [][]arr){
  // ArrayList<ArrayList<Integer>> overinter=new ArrayList<>();
  // int listindexer=0;
  // for(int i=0;i<arr.length;i++){
  // ArrayList<Integer> prerow=new ArrayList<>();
  // if(prerow.get(1)>=arr[i][0]){
  // if()
  // }

  // }
  // }
  public static void mergeTwoSortedArraysWithoutExtraSpace(long[] a, long[] b) {
    // Write your code here.
    int i = a.length - 1, j = 0;
    while (i >= 0 && j < b.length) {
      if (a[i] > b[j]) {
        long temp = b[j];
        b[j] = a[i];
        a[i] = temp;
        i--;
        j++;
      } else {
        break;
      }
    }

    Arrays.sort(a);
    Arrays.sort(b);
  }

  public static boolean rosechecker(int[] a, int mid, int m, int k) {
    int ans = 0;
    int count = 0;
    for (int i = 0; i < a.length; i++) {
      if (mid >= a[i]) {
        count++;
      } else {
        ans += (count / k);
        count = 0;
      }
    }
    ans += (count / k);
    if (ans >= m) {
      return true;
    }
    return false;

  }

  public static int[] getlowest(int[] a) {
    int min = Integer.MAX_VALUE;
    int max = Integer.MIN_VALUE;
    for (int i = 0; i < a.length; i++) {
      max = Math.max(max, a[i]);
      min = Math.min(min, a[i]);
    }
    return new int[] { min, max };
  }

  public static int roseGarden(int[] arr, int m, int k) {
    if (arr.length < (m * k)) {
      return -1;
    }
    int ans = 0;
    int[] minmax = getlowest(arr);
    int low = minmax[0];
    int high = minmax[1];
    while (low <= high) {
      int mid = (low + high) / 2;
      boolean check = rosechecker(arr, mid, m, k);
      if (check) {
        ans = mid;
        high = mid - 1;
      } else {
        low = mid + 1;
      }
    }
    return low;
  }

  public static int checker(int[] a, int mid) {
    int currentlimit = 0;
    for (int i = 0; i < a.length; i++) {
      currentlimit += Math.ceil((double) a[i] / (double) mid);
    }
    return currentlimit;
  }

  public static int smallestDivisor(int arr[], int limit) {
    int low = 1, high = limit;
    int ans = 0;
    while (low <= high) {
      int mid = (low + high) / 2;
      int currentlimit = checker(arr, mid);
      if (currentlimit <= limit) {
        ans = mid;
        high = mid - 1;
      } else {
        low = mid + 1;
      }
    }
    return ans;
  }

  public static int capacitychecker(int[] weigths, int currentcapacity) {
    int days = 1;
    int sum = 0;
    for (int i = 0; i < weigths.length; i++) {
      sum += weigths[i];
      if (sum > currentcapacity) {
        sum = weigths[i];
        days++;
      }
    }
    return days;
  }

  public static int leastWeightCapacity(int[] weights, int d) {
    int high = 0;
    int ans = 0;
    for (int i = 0; i < weights.length; i++) {
      high += weights[i];
    }
    int low = high / d;
    while (low <= high) {
      int mid = (low + high) / 2;
      int currentdays = capacitychecker(weights, mid);
      if (currentdays <= d) {
        ans = mid;
        high = mid - 1;
      } else {
        low = mid + 1;
      }

    }
    return low;
  }

  public static int rowsum(int[][] grid, int i) {
    int sum = 0;
    for (int j = 0; j <= grid[0].length - 1; j++) {
      if (grid[i][j] == 1) {
        sum++;
      } else {
        sum--;
      }
    }
    return sum;
  }

  public static int colsum(int[][] grid, int j) {
    int sum = 0;
    for (int i = 0; i <= grid.length - 1; i++) {
      if (grid[i][j] == 1) {
        sum++;
      } else {
        sum--;
      }
    }
    return sum;
  }

  public static class pair {
    double gap;
    int index;

    pair(double gap, int index) {
      this.gap = gap;
      this.index = index;
    }
  }

  public static double MinimiseMaxDistance(int[] arr, int K) {
    Queue<pair> pq = new PriorityQueue<>((a, b) -> Double.compare(b.gap, a.gap));
    int[] noofgas = new int[arr.length - 1];
    for (int index = 0; index < arr.length - 1; index++) {
      pq.offer(new pair(arr[index + 1] - arr[index], index));
    }
    for (int i = 1; i <= K; i++) {
      pair curr = pq.poll();
      int index = curr.index;
      // double gap=curr.gap;
      noofgas[index]++;
      double actualdis = arr[index + 1] - arr[index];
      double reduce = actualdis / (double) (noofgas[index] + 1);
      pq.offer(new pair(reduce, index));

    }
    return pq.peek().gap;
  }

  // public static void practiseopps(){
  // String s="10111213";
  // Arrays.sort(s,(a,b)->{a.compare(b);});
  // system.out.println(s);
  // }

  // public static int minDifficulty(int[] jobsd, int d) {
  // if(jobs.length>d){
  // return -1;
  // }
  // // return -1; havent done memoization recursion
  // }

  public static double medianbrute(int[] right, int[] left) {
    int[] arr = new int[left.length + right.length];
    int k = 0, l = 0, r = 0;
    while (r < right.length && l < left.length) {
      if (right[r] < left[l]) {
        arr[k++] = right[r++];
      } else {
        arr[k++] = left[l++];
      }
    }
    while (r < right.length) {
      arr[k++] = right[r++];
    }
    while (l < left.length) {
      arr[k++] = left[l++];
    }
    if (arr.length % 2 == 0) {
      return (double) (arr[(int) (arr.length / 2)] + arr[1 + (int) (arr.length / 2)]) / 2;
    }
    return (double) arr[(int) (arr.length / 2)];
  }

   public static double median(int[] a, int[] b) {
        int asize=a.length;
        int bsize=b.length;
        int totalleft=(asize+bsize+1)/2;
        if(asize>bsize) return median(b,a);
        int low=0;
        int high=asize;
        while(low<=high){
          int mid=(low+high)/2;
          int moreleft=totalleft-mid;
          int l1;
          int l2;
          int r1;
          int r2;
          if(mid>0){
            l1=a[mid-1];
          }
          else{
            l1=Integer.MIN_VALUE;
          }
          if(moreleft>0){
            l2=b[moreleft-1];
          }
          else{
            l2=Integer.MIN_VALUE;
          }
          if(mid<asize){
              r1=a[mid];
          }
          else{
              r1=Integer.MAX_VALUE;
          }
          if(moreleft<bsize){
            r2=b[moreleft];
          }
          else{
              r2=Integer.MAX_VALUE;
          }
          if(l1<=r2 && l2<=r1){
            if((asize+bsize)%2==1){
              return (double)Math.max(l1,l2);
            }
            return (double)(Math.max(l1,l2)+Math.min(r1,r2))/2.0;
          }
          if(l1>r2){
            high=mid-1;
          }
          else{
            low=mid+1;
          }
        }
        return 0;
    }


  public static int kthElement(ArrayList<Integer> arr1, ArrayList<Integer> arr2, int k) {
      int count=0;
      int ele=-1;
      int left=0;
      int rigth=0;
      while(left<arr1.size() && rigth<arr2.size()){
        if(count==k){
          return ele;
        }
        if(arr1.get(left)<arr2.get(rigth)){
          ele=arr1.get(left);
          left++;
        }
        else{
          ele=arr2.get(rigth);
          rigth++;
        }
        count++;
      }
      while (left< arr1.size()) {
        if(count==k){
          return ele;
        }
        else{
          ele=arr1.get(left);
          count++;
        }
      }
      while (rigth<arr2.size()) {
         if(count==k){
          return ele;
        }
        else{
          ele=arr2.get(rigth);
          count++;
        }
      }
      return 0;
    
  }

   public static int minSteps(String s, String t) {
        int count=0;
        int i=0;
        int j=0;
        int n=s.length();
        char[] charArray = s.toCharArray();
         char[] charArray1 = t.toCharArray();
        // Sorting the character array
        java.util.Arrays.sort(charArray);
         java.util.Arrays.sort(charArray1);
        // Creating a new string from the sorted character array
        String ss = new String(charArray);
        String tt = new String(charArray1);
        while(i<n && j<n){
          int snum=(int)ss.charAt(i);
          int tnum=(int)tt.charAt(j);
          if(snum==tnum){
            i++;
            j++;
          }
          else if(snum<tnum){
            i++;
            count++;
          }
          else{
            j++;
            count++;
          }
        }
        if(i<n){
          count+=n-i;
        }
        if(j<n){
          count+=n-j;
        }
        return (int)count/2;
    }
    public static int optimalminSteps(String s, String t) {
       int count=0;
    //    int n=s.length();
       int[] scount=new int[26];
       int[] tcount=new int[26];
        for(char ch : s.toCharArray()){
            scount[ch-'a']++;
        }
        for(char ch : t.toCharArray()){
            tcount[ch-'a']++;
        }
        for(int i=0;i<26;i++){
            count+=Math.abs(scount[i]-tcount[i]);
        }
        return count/2;
    }
    
    
    // public boolean closeStrings(String word1, String word2) {
    //     int n=word1.length();
    //     if(n!=word2.length()) return false;
    //     Map<Character,Integer> word1counter=new Treemap<>();
    //     Map<Character,Integer> word2counter=new Treemap<>();
    //     for(int i=0;i<=n-1;i++){
    //       char c=word1.charAt(i);
    //       char e=word2.charAt(i);
    //       if(!word2.contains(i)) return false;
    //       if(word1counter.containsKey(c)){         this apporach is crt but implementation isnt
    //         word1counter.put(c,word1counter.get(c));
    //       }
    //       else{                dont use tree map or map for alphabets because they are only 26
                                  // using arrays it will be useful and easy and effecient
    //         word1counter.put(c,1);
    //       }
    //       if(word2counter.containsKey(e)){
    //         word2counter.put(e,word1counter.get(e));
    //       }
    //       else{
    //         word2counter.put(e,1);
    //       }

    //     }
    //     return false;
    // }


    public static boolean closeStrings(String word1, String word2) {
      if(word1.length()!=word2.length()) return false;
      int n=word1.length();
      int[] freq1=new int[26];
      int[] freq2=new int[26];
      for(int i=0;i<=n-1;i++){
        char c1=word1.charAt(i);
        char c2=word2.charAt(i);
        freq1[c1-'a']++;
        freq2[c2-'a']++;
      }
      for(int i=0;i<=n-1;i++){
            if((freq1[i]==0 && freq2[i]!=0 ) || (freq2[i]==0 && freq1[i]!=0)){
              return false;
            }
      }
      Arrays.sort(freq1);
      Arrays.sort(freq2);
      for(int i=0;i<=n-1;i++){
        if(freq1[i]!=freq2[i]){
          return false;
        }
      }
      return true;


    }

    public static String reverseWordswithotsplit(String s) {
        StringBuilder rev=new StringBuilder();
        for(int i=s.length()-1;i>=0;i--){
          if(s.charAt(i)==' '){
              if(rev.charAt(rev.length()-1)==' '){
                continue;
              }
              else{
                rev.append(' ');
              }
          }
          else{
            rev.append(s.charAt(i));
          }
        }
        String reve=rev.toString();
        Stack<Character> orderch=new Stack<>();
        StringBuilder ans=new StringBuilder();
        for(int i=0;i<reve.length();i++){
            if(reve.charAt(i)==' '){
              ans.append(' ');
              while(!orderch.isEmpty()){
                ans.append(orderch.pop());
              }
            }
            else{
              orderch.push(reve.charAt(i));
            }
        }
        if(!orderch.isEmpty()){
        ans.append(' ');
        while(!orderch.isEmpty()){
          ans.append(orderch.pop());
        }
        }
        return ans.toString().substring(1,ans.length());
    }
    public static String reverseWords(String s) {
        String[] each=s.trim().split("\\s+");
        StringBuilder rev=new StringBuilder();
        for(int i=each.length-1;i>=0;i--){
          rev.append(each[i]);
          rev.append(" ");
        }
        return rev.toString();
    }

    public static String largestOddNumber(String num) {
        for(int i=num.length()-1;i>=0;i--){
          int curpos=Integer.parseInt(""+num.charAt(i));
          if(curpos%2==1){
            return num.substring(0,i+1);
          }
        }
        return "";
    }

    public static String longestCommonPrefix(String[] strs) {
      int minind=-1;
      int min=Integer.MIN_VALUE;
      for(int i=0;i<strs.length;i++){
          if(min<strs[i].length()){
            min=strs[i].length();
            minind=i;
          }
      }
      String minStr=strs[minind];
      for(int i=minStr.length()-1;i>=0;i--){
        String curpre=minStr.substring(0,i+1);
        for(int k=0;k<strs.length;k++){
          if(!strs[k].startsWith(curpre)){
            break;
          }
          else if(strs[k].startsWith(curpre) && k==strs.length-1){
            return curpre;
          }
        }
      }
      return "";
    }


     public static boolean isIsomorphic(String s, String t) {
        int[] smap=new int[200];
        int[] tmap=new int[200];
        for(int i=0;i<s.length();i++){
          if(smap[s.charAt(i)]!=tmap[t.charAt(i)]){
            return false;
          }
        //   int id = (int) (Math.random() * 200) + 1;
          smap[s.charAt(i)]=i+1;
          tmap[t.charAt(i)]=i+1;
        }
        return true;
    }

    public static boolean rotateString(String s, String goal) {
        if(s.length()!=goal.length()) return false;
        char first=s.charAt(0);
        int index=-1;
        for(int i=0;i<s.length();i++){
            if(first==goal.charAt(i)){
              index=i;
              break;
            }
        }
        System.out.print(goal.substring(0,index));
        if(s.contains(goal.substring(0,index)) && s.contains(goal.substring(index,goal.length()))) return true;
        return false;//deonst pass all test case;
        
    }
    static class pairfe{
          int c;
          int freq;
         public pairfe(int c,int freq){
              this.c=c;
              this.freq=freq;
          }
    }
    public static String frequencySort(String s) {
          pairfe[] list=new pairfe[200];
          // System.out.print(list[0]);
          for(int i=0;i<s.length();i++){
              int inde=s.charAt(i);
              if(list[inde]==null){
                list[inde]=new pairfe(inde,1);
              }
              else{
                list[inde].freq++;
              }
          }
          Arrays.sort(list, (a, b) -> {
    if (a == null && b == null) {
        return 0;
    } else if (a == null) {
        return Integer.compare(b.freq, 0);
    } else if (b == null) {
        return Integer.compare(0, a.freq);
    } else {
        return Integer.compare(b.freq, a.freq);
    }
});
StringBuilder str=new StringBuilder();
    for(int i=0;;i++){
      if(list[i]==null){
        break;
      }
      else{
        int co=list[i].freq;
        while(co>0){
            str.append((char)list[i].c);
            co--;
        }
      }
    }
          
          
          return str.toString();
    }

    public static int maxDepth(String s) {
        int maxdep=0;
        int count=0;
        for(int i=0;i<s.length();i++){
          char c=s.charAt(i);
          if(c=='(') count++;
          else if(c==')') count--;
          maxdep=Math.max(count,maxdep);
        }
        return maxdep;
    }    
    
    // public int myAtoi(String s) {
        // StringBuilder number=new StringBuilder();
        // boolean neg=false;
        // if(s.contains("-")) neg=true; 
        // for(int i=0;i<s.length();i++){
          // char cur=s.charAt(i);
          // if(Character.isDigit(cur){
            // 
          // }
          // else{
            // break;
          // }
        // }
        // if(Integer.parseInt(number.toString())>Integer.MAX_VALUE && neg) return Integer.MIN_VALUE;
        // if(In)
    // }
    public static int countSubStrings(String str, int k) {
    HashMap<Character, Integer> hs = new HashMap<>();
    int tsub = 0;
    int i = 0;
    int j = 0;
    while (j <= str.length() - 1) {
        char cur = str.charAt(j);
        if (hs.size() == k) {
            tsub++;

            char ati = str.charAt(i);
            if (hs.get(ati) > 1) {
                hs.put(ati, hs.get(ati) - 1);
            } else {
                hs.remove(ati);
            }
            i++;
        } else if (hs.size() < k) {
            if (hs.containsKey(cur)) {
                hs.put(cur, hs.get(cur) + 1);
            } else {
                hs.put(cur, 1);
            }
            j++;
        }
    }
    // Check the remaining substrings from the last index i to the end
    while (i < str.length() && hs.size() == k) {
        tsub++;
        char ati = str.charAt(i);
        if (hs.get(ati) > 1) {
            hs.put(ati, hs.get(ati) - 1);
        } else {
            hs.remove(ati);
        }
        i++;
    }
    return tsub+1;
  }

  public static String longestPalindromeeror(String s) {//time complexity O(n)
      int i=0;//doesnt pass all the testcases only 72 testcases passed
      int j=s.length()-1;
      if(s.length()==1) return s;
      if(s.length()==2  && s.charAt(i)!=s.charAt(j)) return s.charAt(i)+"";
      // "babad"     "cbbd"
      while(i<=j){
        char val1=s.charAt(i);
        char val2=s.charAt(j);
        if(val2!=val1){
          // System.out.println(val1 +"" + val2);
          System.out.println(s.substring(i,j+1));
          if(s.substring(i+1,j).contains(""+val1)){
             j--;
          }
          else if(s.substring(i,j).contains(""+val2) ) {
            i++;
          }
          else{
            i++;
            j--;
          }
        }
        else{
          int ind1=i;
          int ind2=j;
          while(ind1<=ind2){
            char val3=s.charAt(ind1);
            char val4=s.charAt(ind2);
            if(val3!=val4) {
              i=ind1;
              j=ind2;
              break;
            }
            else if(ind1==ind2 || ind1+1==ind2){
               return s.substring(i,j+1);
            }
            System.out.println(s.substring(ind1,ind2+1));
            ind1++;
            ind2--;
          }
        }
      }
      return ""+s.charAt(0);                    
  }
  public static int evaluate(int num1,int num2,char cur){
      switch(cur){
        case '+' : return num1+num2;
        case '-' : return num1-num2;
        case '/' : return (int)(num1/num2);
        case '*' : return num1*num2;
      }
      return 0;
  }

  public static int evalRPN(String[] tokens) {
        Stack<Integer> st=new Stack<>();
        for(int i=0;i<tokens.length;i++){
          char cur=tokens[i].charAt(tokens[i].length()-1);
          if(Character.isDigit(cur)){
            st.push(Integer.parseInt(tokens[i]));
            // System.out.println(tokens[i]);
          }
          else{
            int num2=st.pop();
            int num1=st.pop();
            st.push(evaluate(num1,num2,cur));
            // System.out.println(st.peek());
          }
        } 
        return st.peek();     
  }

  public static String longestPalindrome(String s){
      String res="";
      for(int i=0;i<s.length;i++){
       }
  
  }


  public static void main(String[] args) {
        int[] arr={1,2,3,4,5,6,7,8,9,10};
        // arr=selection(arr);
        /*  for (int i = 0; i < arr.length; i++) {
            System.out.println(arr[i]);
        }*/
        // String k=removeConsecutiveDuplicates("aabbcb");`
            // List<List<Integer>> l=generate(11);
            // for (List<Integer> row : l) {
            // for (int element : row) {
            //     System.out.print(element + " ");
            // }
            // System.out.println(); // Move to the next line for a new row
            // }
            // if(false){
            //     System.out.println("if");
            // }
            // else if(true){
            //     System.out.println("else if");
            // }
            // else if(true){
            //     System.out.println("else 2if");
            // }
            // else{
            //     System.out.println("hello");
            // }

            // int[][] e=;
            // int[] r=spiralmatrix(e);
            // triplet(4,arr);
            // System.out.println(leastWeightCapacity(arr,5));
        // int[][] di={{0,1,1},{1,0,1},{0,0,1}};
        // int[][] diff=onesMinusZeros(di);
        // for(int i=0;i<=diff.length-1;i++){
            // for(int j=0;j<diff[0].length;j++){
                // System.out.println(diff[i][j]);
            // }
            // int[] ar={1,2,1,0,0,0,2};
        // }
            // sort012(ar);
            // int[] a = {1,10};
        // int k = 10;
        // double ans = MinimiseMaxDistance(a, k);
        // System.out.println("The answer is: " + ans);
        // int k=1000;
        // int l=1000;
        // System.out.println(k==l);
        // practiseopps();
        // int[] a = {1, 4, 7, 10, 12};
        // int[] b = {2, 3, 6, 15};
        //   System.out.println(median(a,b));
      // System.out.println(8);
      // ArrayList<Integer> l=new ArrayList<>();
      // ArrayList<Integer> l2=new ArrayList<>();
      // l.add(1);
      // l.add(3);
      // l.add(30);
      // l2.add(12);
      // System.out.print(kthElement(l,l2,4));
      // System.out.print(minSteps("leetcode","practice"));
      // System.out.print(closeScd trings("abbccc","cabbba"));
        // System.out.println(reverseWords("hello world to be     happy"));
        // String s="   dsd m    kdm    ";
        // String[] sds=s.trim().split("\\s+");
        // System.out.println(largestOddNumber("52"));
        // String[] strs = {"flower","flow","flight"};
        // System.out.println(longestCommonPrefix(strs));
        // System.out.println(isIsomorphic("fok","hdd"));
        // String s="nds";
        // String goal="sd";
        // System.out.println(rotateString("bbbacddceeb","ceebbbbacdd"));
        // System.out.println(frequencySort("rttetrrreer"));
          // System.out.println(maxDepth("(1)+((2))+(((3)))"));
          // System.out.print(countSubStrings("aacfssa",3));
          // String[] tokens={"10","6","9","3","+","-11","*","/","*","17","+","5","+"};
          // System.out.println(evalRPN(tokens));
          System.out.println(longestPalindrome("aacabdkacaa"));
      }

}
