/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package AoUtils;

import format.position.ChrPos;
import utils.IOUtils;
import utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.text.SimpleDateFormat;
import java.util.*;

/**
 *
 * @author Aoyue
 */
public class Mode {

    /*==================================== 测试用 =============================================*/
    public Mode() {

    }

    /**
     *
     * @param infileDirS
     * @param outfileDirS
     */
    public void vcfParallel(String infileDirS, String outfileDirS, String extractRatio) {

        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        Collections.sort(fsList);
        fsList.parallelStream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                String outfileS = null;
                BufferedReader br = null;
                if (infileS.endsWith(".vcf")) {
                    br = IOUtils.getTextReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().split(".vcf")[0] + "_subset.vcf.gz").getAbsolutePath();
                } else if (infileS.endsWith(".vcf.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().split(".vcf.gz")[0] + "_subset.vcf.gz").getAbsolutePath();
                }
                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                String temp = null;
                int cnttotal = 0;
                int cntsubset = 0;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) {
                        bw.write(temp);
                        bw.newLine();
                    } else {
                        cnttotal++;
                        double r = Math.random();
                        double ratio = Double.parseDouble(extractRatio);
                        if (r > ratio) {
                            continue; //返回带正号的 double 值，该值大于等于 0.0 且小于 1.0。返回值是一个伪随机选择的数，在该范围内（近似）均匀分布
                        }
                        List<String> l = PStringUtils.fastSplit(temp);
                        if (l.get(3).contains(",")) {
                            continue; // 第3列是alt的信息，若有2个等位基因，则去除这一行
                        }
                        bw.write(temp);
                        bw.newLine();
                        cntsubset++;
                    }
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(f.getName() + "\twith " + cnttotal + " bp has a subset of\t" + cntsubset + "\tis completed at " + outfileS);
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    //本模板使用多线程流
    public void txtParallel() {
        String infileDirS = "";
        String outfileDirS = "";
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        Collections.sort(fsList);
        fsList.parallelStream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                String outfileS = null;
                BufferedReader br = null;
                if (infileS.endsWith(".txt")) {
                    br = IOUtils.getTextReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().split(".txt")[0] + "_subset.txt.gz").getAbsolutePath();
                } else if (infileS.endsWith(".txt.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().split(".txt.gz")[0] + "_subset.txt.gz").getAbsolutePath();
                }
                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    List<String> l = PStringUtils.fastSplit(temp);

                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(f.getName() + "\tis completed at " + outfileS);
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    public void mode() {
        String infileDirS = "";
        String outfileDirS = "";

        //  new CountSites().mergesubsetVCF(args[0], args[1]);
        // new CountSites().mergesubsetVCF(args[0], args[1], args[2]);
/*==================================== 测试用 =============================================*/
        try {
            String infileS = "";
            String outfileS = "";
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = null;
            int cnt = 0;
            List<String> l = new ArrayList<>();
            while ((temp = br.readLine()) != null) {
                StringBuilder sb = new StringBuilder();
                if (temp.startsWith("#")) {
                    continue;
                }
                cnt++;
            }
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        /*==================================== 测试用 =============================================*/
        try {
            String infileS = "";
            String outfileS = "";
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);

            StringBuilder sb = new StringBuilder();
            sb.append("##fileformat=VCFv4.1\n");
            SimpleDateFormat sdf = new SimpleDateFormat("MM/dd/yyyy HH:mm:ss.SSS");
            Date dt = new Date();
            String S = sdf.format(dt);
            sb.append("##fileDate=").append(S.split(" ")[0]).append("\n");
            //sb.append("##reference=").append(referenceFileS).append("\n");
            sb.append("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"").append("Total depth").append("\">\n");
            sb.append("##INFO=<ID=AD,Number=2+,Type=Integer,Description=\"").append("Total allelelic depths in order listed starting with REF").append("\">\n");
            sb.append("##INFO=<ID=NZ,Number=1,Type=Integer,Description=\"").append("Number of individuals with alleles present").append("\">\n");
            sb.append("##INFO=<ID=AP,Number=2+,Type=Integer,Description=\"").append("Number of individuals in which an allele is present").append("\">\n");
            sb.append("##INFO=<ID=PV,Number=1+,Type=Float,Description=\"").append("Segreagation test P-Value of alternative alleles").append("\">\n");
            sb.append("##INFO=<ID=DI,Number=2,Type=Integer,Description=\"").append("Number of deletion and insertion type").append("\">\n");
            sb.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"").append("Genotype").append("\">\n");
            sb.append("##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"").append("Allelic depths for the reference and alternate alleles in the order listed").append("\">\n");
            sb.append("##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"").append("Genotype likelihoods for 0/0, 0/1, 1/1, or 0/0, 0/1, 0/2, 1/1, 1/2, 2/2 if 2 alt alleles").append("\">\n");
            //return sb.toString();

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        //////////////////////
        try {
            String infileS = "";
            String outfileS = "";
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                cnt++;

            }
            br.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

//////////////////////////////////////////////////////////////////
        List<ChrPos> l = new ArrayList();
        String chr = null;
        String pos = null;
        l.add(new ChrPos(Short.valueOf(chr), Integer.valueOf(pos)));

    }

    public void testifD() {
        boolean ifd = false;
        String chr = "chr003.vmap2...".substring(3, 6);
        //根据染色体号进行AB还是D的判断
        String[] db = {"5", "6", "11", "12", "17", "18", "23", "24", "29", "30", "35", "36", "41", "42"};
        Arrays.sort(db);
        if (Arrays.binarySearch(db, chr) > -1) { //说明是属于D的
            ifd = true;
        }
        try {
            String infileS = "";
            String outfileS = "";
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;

            if (ifd == false) {
                bw.write(temp + "\tAncestral\tDaf\tDaf_ABD\tDaf_AB");
                bw.newLine();
            } else if (ifd == true) {
                bw.write(temp + "\tAncestral\tDaf\tDaf_ABD\tDaf_D");
                bw.newLine();
            }

            while ((temp = br.readLine()) != null) {
                cnt++;

            }
            br.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }
    
    /**
     * 在输出目录中建立相同的子目录
     */
    public void mksamedirs() {
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/001_ori";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/002_changeChrPos";
        File[] fs = new File(infileDirS).listFiles();
        for(int i = 0; i < fs.length; i++){
            if(fs[i].isHidden())
                fs[i].delete();
        }
        fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            System.out.println(fs[i]);
            new File(outfileDirS, fs[i].getName()).mkdirs();
        }
        
        for(int i = 0; i < fs.length; i++){
            new CountSites().mergefileandChangeChrPos_chr1and2(fs[i].getAbsolutePath(),new File(outfileDirS, fs[i].getName()).getAbsolutePath());
        }
    }

}
