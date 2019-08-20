/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package AoUtils;

import format.table.RowTable;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.Arrays;
import java.util.List;
import utils.IOUtils;
import utils.PArrayUtils;
import utils.PStringUtils;

/**
 *
 * @author Aoyue
 */
public class Bin {

    public Bin() {

    }
    
    /**
     * 根据合并后的1D 2D 3D 4D 5D 6D 7D文件，分别建立BinTable,最终依旧在一个文件中输出。
     * @param infileS
     * @param outfileDirS
     * @param binSize 
     */
    public void mkBinTable(String infileS, String outfileDirS, int binSize){
        
    }

    
    /**
     * make bin table, 一个文件制作一个Binfile
     * @param infileDirS
     * @param outfileDirS
     * @param binSize 
     */
    public void mkBintable(String infileDirS, String outfileDirS, int binSize) {
        //String 
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        fsList.stream().forEach(f -> {
            String infileS = f.getAbsolutePath();
            String outfileS = null;
            BufferedReader br = null;
            if (infileS.endsWith(".txt")) {
                br = IOUtils.getTextReader(infileS);
                outfileS = new File(outfileDirS, f.getName().replaceFirst(".txt", "." + binSize/1000000 + "M" + ".binTable.txt")).getAbsolutePath();
            } else if (infileS.endsWith(".txt.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
                outfileS = new File(outfileDirS, f.getName().replaceFirst(".txt.gz", "." + binSize/1000000 + "M" + ".binTable.txt")).getAbsolutePath();
            }

            RowTable<String> t = new RowTable(infileS);
            //染色体的长度是最后一行pos的位置
            int chrlength = Integer.valueOf(t.getCell(t.getRowNumber() - 1, 1));
            String Chr = t.getCell(0, 0);

            int[][] bound = PArrayUtils.getSubsetsIndicesBySubsetSize(chrlength, binSize);
            int count[] = new int[bound.length];
            int[] bounds = new int[bound.length];
            for(int i=0; i< bound.length; i++){
                bounds[i] = bound[i][0];
            }
            for (int i = 0; i < t.getRowNumber(); i++) {
                int index = Arrays.binarySearch(bounds, Integer.valueOf(t.getCell(i, 1)));
                if (index < 0) {
                    index = -index - 2;
                }
                count[index]++;
            }

            try {
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                bw.write("CHROM\tBIN_START\tSNP_COUNT\tVARIANTS.KB");
                bw.newLine();
                for (int i = 0; i < bound.length; i++) {
                    StringBuilder sb = new StringBuilder();
                    double variant = (double)count[i] / (double)1000;
                    sb.append(Chr).append("\t").append(bound[i][0]).append("\t").append(count[i]).append("\t").append(variant);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();

            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
    }

    /**
     * extract pos info from vcf file. eg: vcf ---- Chr Pos
     *
     * @param infileDirS
     * @param outfileDirS
     */
    public void mkHapPos(String infileDirS, String outfileDirS) {
        File[] fs = new File(infileDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        fsList.stream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                String outfileS = null;
                BufferedReader br = null;
                if (infileS.endsWith(".vcf")) {
                    br = IOUtils.getTextReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().replaceFirst(".vcf", ".pos.txt.gz")).getAbsolutePath();
                } else if (infileS.endsWith(".vcf.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().replaceFirst(".vcf.gz", ".pos.txt.gz")).getAbsolutePath();
                }

                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                String temp = null;
                int cnt = 0;
                bw.write("Chr\tPos\n");
                List<String> l = null;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) {
                        continue;
                    }
                    temp = temp.substring(0, 40); //肯定够                 
                    l = PStringUtils.fastSplit(temp);
                    StringBuilder sb = new StringBuilder();
                    sb = new StringBuilder(l.get(0));
                    sb.append("\t").append(l.get(1));
                    bw.write(sb.toString());
                    bw.newLine();
                    if (cnt % 1000 == 0) {
                        System.out.println("Output " + String.valueOf(cnt) + " SNPs");
                    }
                    cnt++;
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(String.valueOf(cnt) + " SNPs output from " + f.getAbsolutePath());

            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
    }

}
