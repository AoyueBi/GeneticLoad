/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package WheatGeneticLoad;

import pgl.format.table.RowTable;
import gnu.trove.list.array.TIntArrayList;
import pgl.utils.IOUtils;
import pgl.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import static java.util.stream.Collectors.toList;

/**
 *
 * @author Aoyue
 */
public class TreePreparation {

    public TreePreparation() {
//        this.getintersection();
        this.getsubsetFasta();

    }







    /**
     * 改变方针策略：先进行随机抽取barleyde 文件的pos和ATGC，然后在vmap2里边进行寻找
     *
     */

    /**
     * 根据posList库，进行barley的抽样和VCF的抽样
     *
     */
    public void getsubsetFasta(){
        String infileDirS = "/Users/Aoyue/Documents/test/Vmap2";
        String posDirS = "/Users/Aoyue/Documents/test/out1";
        String ancesDirS = "/Users/Aoyue/Documents/test/AncesbyChr"; //ancestral allele的路径
        String outfileDirS = "/Users/Aoyue/Documents/test/out2";
        String vcfoutDirS = "/Users/Aoyue/Documents/test/out3";
        double ratio = 0.5;
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
                //处理ancestral数据库的文件路径
                String infileS = f.getAbsolutePath();
                String outfileS = null;
                String posS = null;
                String ancS = null;
                BufferedReader br = null;
                if (infileS.endsWith(".vcf")) {
                    outfileS = new File(outfileDirS, f.getName().replaceFirst("_vmap2.1.vcf", "_barley.fasta.txt.gz")).getAbsolutePath();
                    posS = new File(posDirS, f.getName().replaceFirst("_vmap2.1.vcf", "_barley.intersection.txt.gz")).getAbsolutePath();
                    ancS = new File(ancesDirS, f.getName().replaceFirst("_vmap2.1.vcf", "_ancestral.txt.gz")).getAbsolutePath();
                } else if (infileS.endsWith(".vcf.gz")) {
                    outfileS = new File(outfileDirS, f.getName().replaceFirst("_vmap2.1.vcf.gz", "_barley.fasta.txt.gz")).getAbsolutePath();
                    posS = new File(posDirS, f.getName().replaceFirst("_vmap2.1.vcf.gz", "_barley.intersection.txt.gz")).getAbsolutePath();
                }


                //先建立的posList1 进行抽样
                TIntArrayList posList = new TIntArrayList();
                br = IOUtils.getTextGzipReader(posS);
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    double r = Math.random();
                    if (r > ratio) {
                        continue; //返回带正号的 double 值，该值大于等于 0.0 且小于 1.0。返回值是一个伪随机选择的数，在该范围内（近似）均匀分布
                    }
                    posList.add(Integer.parseInt(temp));
                }
                br.close();

                //输出文件是barley的fasta格式
                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS); //输出文件是压缩格式的
                bw.write(">Barley");
                bw.newLine();
                br = IOUtils.getTextGzipReader(ancS);
                br.readLine();
                List<String> l = new ArrayList<>();
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    int pos = Integer.parseInt(l.get(1));
                    String barley = l.get(5);
                    int index = posList.binarySearch(pos);
                    if(index > -1){
                        bw.write(barley);
                    }
                }
                bw.newLine();
                br.close();
                bw.flush();bw.close();

                br=IOUtils.getTextReader(infileS);
                String vcfS = new File(vcfoutDirS,f.getName().replaceFirst("_vmap2.1.vcf", "_vmap2.1_subset_.vcf.gz")).getAbsolutePath();
                bw = IOUtils.getTextGzipWriter(vcfS);
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) {
                        bw.write(temp);
                        bw.newLine();
                    } else {
                        l = PStringUtils.fastSplit(temp.substring(0,25));
                        int pos = Integer.parseInt(l.get(1));
                        int index = posList.binarySearch(pos);
                        if(index > -1){
                            bw.write(String.valueOf(temp));
                            bw.newLine();
                        }
                    }
                }
                bw.flush();
                bw.close();
                br.close();

            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    /**
     *致命性缺点：程序已放弃，吐核
     *
     */
//    画树步骤：
//            1.获取Chr	Pos	Barley 子集；
//            2.写一个方法，根据posList获取vcf文件子集，返回的是VCF文件；
//            3.写一个方法，将第一步的结果转化成fasta格式。
    // outfileS = new File(outfileDirS, f.getName().replaceFirst("_ancestral.txt","_vmap2.1.vcf")).getAbsolutePath();
    public void getintersection() {

//        String  ancesDirS = "/Users/Aoyue/Documents/test/AncesbyChr"; //ancestral allele的路径
//        String outfileDirS = "/Users/Aoyue/Documents/test/out1"; //输出chr pos barley文件的路径
//        String infileDirS = "/Users/Aoyue/Documents/test/Vmap2";
        String ancesDirS = "/data4/home/aoyue/vmap2/feilu/003_annotation/002_ancestral/byChr"; //ancestral allele的路径
        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/016_tree/001_barleyequence/001_intersection"; //输出chr pos barley文件的路径
        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/011_VMapII";

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
                //处理ancestral数据库的文件路径
                String infileS = f.getAbsolutePath();
                String outfileS = null;
                String ancS = null;
                BufferedReader br = null;
                if (infileS.endsWith(".vcf")) {
                    br = IOUtils.getTextReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().replaceFirst("_vmap2.1.vcf", "_barley.intersection.txt.gz")).getAbsolutePath();
                    ancS = new File(ancesDirS, f.getName().replaceFirst("_vmap2.1.vcf", "_ancestral.txt.gz")).getAbsolutePath();
                } else if (infileS.endsWith(".vcf.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().replaceFirst("_vmap2.1.vcf.gz", "_barley.intersection.txt.gz")).getAbsolutePath();
                    ancS = new File(ancesDirS, f.getName().replaceFirst("_vmap2.1.vcf.gz", "_ancestral.txt.gz")).getAbsolutePath();
                }
                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS); //输出文件是压缩格式的

                //先建立barley的posList1
                List<Integer> posList1 = new ArrayList<>();
                List<Integer> posList2 = new ArrayList<>();
                List<Integer> posList3 = new ArrayList<>();
                RowTable<String> t = new RowTable<>(ancS);
                for (int i = 0; i < t.getRowNumber(); i++) {
                    int pos = Integer.parseInt(t.getCell(i, 1));
                    String ref = t.getCell(i, 3);
                    String barley = t.getCell(i, 5);
                    if (barley.equals("-")) {
                        continue;
                    }
                    if (ref.startsWith("N")) {
                        continue;
                    }
                    posList1.add(pos);
                }
                Collections.sort(posList1);
                System.out.println(posList1.size() + " bp in barley " + f.getName().substring(3,6));

                String temp = null;
                List<String> l = new ArrayList<>();

                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) {

                    } else { //
                        temp = temp.substring(0, 20);
                        l = PStringUtils.fastSplit(temp);
                        int pos = Integer.parseInt(l.get(1));
                        posList2.add(pos);
                    }//
                }
                br.close();
                Collections.sort(posList1);
                System.out.println(posList2.size() + " bp in vmap2 " + f.getName().substring(3, 6));

                posList3 = posList1.stream().filter(item -> posList2.contains(item)).collect(toList());
                Collections.sort(posList3);
                System.out.println(posList3.size() + " bp intersection in file " + f.getName());

                for (int i = 0; i < posList3.size(); i++) {
                    bw.write(String.valueOf(posList3.get(i)));
                    bw.newLine();
                }

                bw.flush();
                bw.close();

//                System.out.println(f.getName() + "\twith " + cnttotal + " bp has a subset of\t" + cntsubset + "\tbiallelic SNPs is completed at " + outfileS);
            } catch (Exception e) {
                e.printStackTrace();
            }
        });

    }

}
