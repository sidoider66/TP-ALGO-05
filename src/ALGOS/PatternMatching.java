package ALGOS;

import java.io.*;
import java.nio.file.*;
import java.util.*;
import javax.swing.*;
import javax.swing.table.DefaultTableModel;

public class PatternMatching {

    public static void main(String[] args) {
        try {
            String filePath = "C:\\Users\\THINKPAD\\OneDrive\\Bureau\\input.txt";
            String text = new String(Files.readAllBytes(Paths.get(filePath)));
            String pattern = "motif";

            Map<String, Long> executionTimes = new LinkedHashMap<>();
            Map<String, Map<String, Object>> detailedResults = new LinkedHashMap<>();

            long startTime = System.nanoTime();
            List<Integer> naiveResults = naiveSearch(text, pattern);
            long elapsedTime = System.nanoTime() - startTime;
            executionTimes.put("Naïve", elapsedTime);
            detailedResults.put("Naïve", Map.of("Indices", naiveResults));

            startTime = System.nanoTime();
            List<Integer> rkResults = rabinKarpSearch(text, pattern);
            elapsedTime = System.nanoTime() - startTime;
            executionTimes.put("Rabin-Karp", elapsedTime);
            detailedResults.put("Rabin-Karp", Map.of("Indices", rkResults));

            startTime = System.nanoTime();
            int[] lps = new int[pattern.length()];
            List<Integer> kmpResults = kmpSearch(text, pattern, lps);
            elapsedTime = System.nanoTime() - startTime;
            executionTimes.put("KMP", elapsedTime);
            detailedResults.put("KMP", Map.of("Indices", kmpResults, "LPS Table", lps));

            startTime = System.nanoTime();
            BMResult bmResult = boyerMooreSearch(text, pattern);
            elapsedTime = System.nanoTime() - startTime;
            executionTimes.put("Boyer-Moore", elapsedTime);
            detailedResults.put("Boyer-Moore", Map.of(
                "Indices", bmResult.indices,
                "Bad Char Table", bmResult.badChar,
                "Suffix Table", bmResult.suff,
                "Good Suffix Table", bmResult.goodSuffix
            ));

            showDetailedResultsInTables(executionTimes, detailedResults);

        } catch (IOException e) {
            System.err.println("Erreur lors de la lecture du fichier : " + e.getMessage());
        }
    }

    public static List<Integer> naiveSearch(String text, String pattern) {
        List<Integer> indices = new ArrayList<>();
        int n = text.length();
        int m = pattern.length();

        for (int i = 0; i <= n - m; i++) {
            int j;
            for (j = 0; j < m; j++) {
                if (text.charAt(i + j) != pattern.charAt(j)) {
                    break;
                }
            }
            if (j == m) {
                indices.add(i);
            }
        }
        return indices;
    }

    public static List<Integer> rabinKarpSearch(String text, String pattern) {
        List<Integer> indices = new ArrayList<>();
        int n = text.length();
        int m = pattern.length();
        int prime = 101;
        int hashPattern = 0, hashText = 0;
        int h = 1;

        for (int i = 0; i < m - 1; i++) {
            h = (h * 256) % prime;
        }

        for (int i = 0; i < m; i++) {
            hashPattern = (256 * hashPattern + pattern.charAt(i)) % prime;
            hashText = (256 * hashText + text.charAt(i)) % prime;
        }

        for (int i = 0; i <= n - m; i++) {
            if (hashPattern == hashText) {
                int j;
                for (j = 0; j < m; j++) {
                    if (text.charAt(i + j) != pattern.charAt(j)) {
                        break;
                    }
                }
                if (j == m) {
                    indices.add(i);
                }
            }
            if (i < n - m) {
                hashText = (256 * (hashText - text.charAt(i) * h) + text.charAt(i + m)) % prime;
                if (hashText < 0) {
                    hashText += prime;
                }
            }
        }
        return indices;
    }

    public static List<Integer> kmpSearch(String text, String pattern, int[] lps) {
        computeLPS(pattern, lps);
        List<Integer> indices = new ArrayList<>();
        int n = text.length();
        int m = pattern.length();

        int i = 0, j = 0;
        while (i < n) {
            if (pattern.charAt(j) == text.charAt(i)) {
                i++;
                j++;
            }
            if (j == m) {
                indices.add(i - j);
                j = lps[j - 1];
            } else if (i < n && pattern.charAt(j) != text.charAt(i)) {
                if (j != 0) {
                    j = lps[j - 1];
                } else {
                    i++;
                }
            }
        }
        return indices;
    }

    private static void computeLPS(String pattern, int[] lps) {
        int m = pattern.length();
        int length = 0, i = 1;

        while (i < m) {
            if (pattern.charAt(i) == pattern.charAt(length)) {
                length++;
                lps[i] = length;
                i++;
            } else {
                if (length != 0) {
                    length = lps[length - 1];
                } else {
                    lps[i] = 0;
                    i++;
                }
            }
        }
    }

    public static BMResult boyerMooreSearch(String text, String pattern) {
        int n = text.length();
        int m = pattern.length();
        int[] badChar = buildBadCharTable(pattern);
        int[] suff = buildSuffixArray(pattern);
        int[] goodSuffix = buildGoodSuffixTable(pattern, suff);

        List<Integer> indices = new ArrayList<>();
        int shift = 0;
        while (shift <= n - m) {
            int j = m - 1;

            while (j >= 0 && pattern.charAt(j) == text.charAt(shift + j)) {
                j--;
            }

            if (j < 0) {
                indices.add(shift);
                shift += (shift + m < n) ? m - badChar[text.charAt(shift + m)] : 1;
            } else {
                shift += Math.max(1, Math.max(j - badChar[text.charAt(shift + j)], goodSuffix[j]));
            }
        }
        return new BMResult(badChar, suff, goodSuffix, indices);
    }

    private static int[] buildBadCharTable(String pattern) {
        int[] badChar = new int[256];
        Arrays.fill(badChar, -1);
        for (int i = 0; i < pattern.length(); i++) {
            badChar[pattern.charAt(i)] = i;
        }
        return badChar;
    }

    private static int[] buildSuffixArray(String pattern) {
        int m = pattern.length();
        int[] suff = new int[m];
        suff[m - 1] = m;
        for (int i = m - 2, g = m - 1, f = m - 1; i >= 0; i--) {
            if (i > g && suff[i + m - 1 - f] < i - g) {
                suff[i] = suff[i + m - 1 - f];
            } else {
                if (i < g) {
                    g = i;
                }
                f = i;
                while (g >= 0 && pattern.charAt(g) == pattern.charAt(g + m - 1 - f)) {
                    g--;
                }
                suff[i] = f - g;
            }
        }
        return suff;
    }

    private static int[] buildGoodSuffixTable(String pattern, int[] suff) {
        int m = pattern.length();
        int[] goodSuffix = new int[m];
        Arrays.fill(goodSuffix, m);
        for (int i = m - 1, j = 0; i >= 0; i--) {
            if (suff[i] == i + 1) {
                for (; j < m - 1 - i; j++) {
                    if (goodSuffix[j] == m) {
                        goodSuffix[j] = m - 1 - i;
                    }
                }
            }
        }
        for (int i = 0; i < m - 1; i++) {
            goodSuffix[m - 1 - suff[i]] = m - 1 - i;
        }
        return goodSuffix;
    }

    private static void showDetailedResultsInTables(Map<String, Long> executionTimes, Map<String, Map<String, Object>> detailedResults) {
        JFrame frame = new JFrame("Pattern Matching Results");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        String[] columnNames = {"Algorithm", "Execution Time (ns)", "Details"};
        DefaultTableModel model = new DefaultTableModel(columnNames, 0);

        for (Map.Entry<String, Long> entry : executionTimes.entrySet()) {
            String algorithm = entry.getKey();
            Long time = entry.getValue();
            Map<String, Object> details = detailedResults.get(algorithm);

            StringBuilder detailsString = new StringBuilder();
            details.forEach((key, value) -> detailsString.append(key).append(": ").append(value.toString()).append("\n"));

            model.addRow(new Object[]{algorithm, time, detailsString.toString()});
        }

        JTable table = new JTable(model);
        JScrollPane scrollPane = new JScrollPane(table);
        frame.add(scrollPane);
        frame.pack();
        frame.setVisible(true);
    }

    static class BMResult {
        int[] badChar;
        int[] suff;
        int[] goodSuffix;
        List<Integer> indices;

        public BMResult(int[] badChar, int[] suff, int[] goodSuffix, List<Integer> indices) {
            this.badChar = badChar;
            this.suff = suff;
            this.goodSuffix = goodSuffix;
            this.indices = indices;
        }
    }
}
