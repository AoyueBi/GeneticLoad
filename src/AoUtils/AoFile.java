/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package AoUtils;

import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.text.SimpleDateFormat;
import java.util.*;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;

/**
 * @author AoyueBi
 *
 */
public class AoFile {
    public AoFile(){
        
    }

    /**
     * 将结果根据亚基因组(Group)分开,写出N个文件
     */
    public static void splitFilebyGroup(String infileS, int groupIndex, String outfileDirS ){
        File f = new File(infileS);
        String[] group = AoFile.getStringArraybySet(infileS,1);
        Arrays.sort(group);
        String[] outfilesS = new String[group.length];
        for (int i = 0; i < group.length; i++) {
            outfilesS[i] = new File(outfileDirS,f.getName().split(".txt")[0] + "_" + group[i] + ".txt").getAbsolutePath();
        }

        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter[] bw = new BufferedWriter[group.length];
            for (int i = 0; i < bw.length; i++) {
                bw[i] = AoFile.writeFile(outfilesS[i]);
            }
            String header = br.readLine();
            for (int i = 0; i < bw.length; i++) {
                bw[i].write(header);
                bw[i].newLine();
            }
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                String query = l.get(groupIndex);
                int index = Arrays.binarySearch(group,query);
                bw[index].write(temp);
                bw[index].newLine();
            }
            br.close();
            for (int i = 0; i < bw.length; i++) {
                bw[i].flush();
                bw[i].close();
            }

            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    /**
     * 根据输入文件，自动创建一个文件，该文件目录与输入文件的父目录相同，目录名字为A_out,文件名字和输入文件名字相同。
     * @param infileS
     * @return
     */
    public static String mkOutfileS(String infileS){
        String outfileS = null;
        String outfileDirS = new File(infileS).getParent(); //获取输入文件的目录
        outfileDirS = new File(outfileDirS).getParent(); //根据输入文件的目录获取上一级父目录；
        outfileDirS = outfileDirS + "/A_out"; //根据父目录路径 创建输出文件目录
        new File(outfileDirS).mkdirs(); //创建输出文件目录
        outfileS = new File(outfileDirS,new File(infileS).getName()).getAbsolutePath();

        return outfileS;
    }


    /**
     * 提取文件的某几列，并且跳过以 skipString 开头的行
     * @param infileS
     * @param skipString
     * @param columns
     * @param outfileS
     */
    public static void extractFileColumn(String infileS, String skipString, int[] columns,String outfileS){

        Arrays.sort(columns);

        try{
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = null;
            String header = br.readLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                if (temp.startsWith(skipString))continue;
                StringBuilder sb = new StringBuilder();
                for (int i = 0; i < l.size(); i++) {
                    int index = Arrays.binarySearch(columns,i);
                    if (index < 0) continue;
                    sb.append(l.get(i)).append("\t");
                }
                sb.deleteCharAt(sb.length()-1);
                bw.write(sb.toString());
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println();

        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 为文本添加分组，根据文件名字
     */
    public static void addGrouptoFile(String infileS,int IndexBegin,int IndexEnd){

//        String infileS="/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/001_depth/002_depthCal/chr1A_vmap2_subset0.001_depth.txt.gz";
        File f = new File(infileS);
        String outfileS = f.getName().split(".txt")[0] + "_addGroup.txt.gz";
        outfileS= new File(f.getParent(),outfileS).getAbsolutePath();

        String group = f.getName().substring(IndexBegin,IndexEnd);
        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = br.readLine();
            bw.write(temp + "\tGroup");
            bw.newLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                bw.write(temp+ "\t" +group);
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }


    public static File filterTxtLines(String infileS, int columIndex, List<String> inlist, String outfileS){
        File out = new File(outfileS);
        Collections.sort(inlist);
        try {
            BufferedReader br = new AoFile().readFile(infileS);
            BufferedWriter bw = new AoFile().writeFile(outfileS);
            String header = br.readLine();
            bw.write(header); bw.newLine();
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            int cntkeep = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                String query = l.get(columIndex);
                int index = Collections.binarySearch(inlist,query);
                if (index < 0) continue;
                cntkeep++;
                bw.write(temp);
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println(cnt + "\tkeep lines " + cntkeep + "\t" + infileS + " is completed at " + outfileS );
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        return out;
    }


    /**
     *  返回满足条件的某些行，即取子集
     *
     * @param infileS
     * @param columIndex
     * @param inlist
     * @param outfileS
     * @return
     */
    public static File filterTxtLines(String infileS, int columIndex, TIntArrayList inlist, String outfileS){
        File out = new File(outfileS);
        inlist.sort();
        try {
            BufferedReader br = new AoFile().readFile(infileS);
            BufferedWriter bw = new AoFile().writeFile(outfileS);
            String header = br.readLine();
            bw.write(header); bw.newLine();
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            int cntkeep = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                int pos = Integer.parseInt(l.get(columIndex));
                int index = inlist.binarySearch(pos);
                if (index < 0) continue;
                cntkeep++;
                bw.write(temp);
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println(cnt + "\tkeep lines " + cntkeep + "\t" + infileS + " is completed at " + outfileS );
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        return out;
    }

    /**
     * 获取表格的列数，不包括表头
     * @param infileS
     * @return
     */
    public static int countFileColumnNumber(String infileS){
        int cntColumn = 0;
        try{
            BufferedReader br = AoFile.readFile(infileS);
            String temp = br.readLine();
            cntColumn = PStringUtils.fastSplit(temp).size();

        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        return cntColumn;
    }
    /**
     * 获取表格的行数，不包括表头
     * @param infileS
     * @return
     */
    public static int countFileRowNumber_withHeader(String infileS){
        int out = 0;
        try{
            BufferedReader br = AoFile.readFile(infileS);
            String temp = br.readLine();
            int cnt=0;
            while((temp=br.readLine()) != null){
                cnt++;
            }
            br.close();
            out = cnt;

        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        return out;
    }

    /**
     * 获取文件数组
     * @param inDirS
     * @return
     */
    public static File[] getFileArrayInDir (String inDirS) {
        File[] fs = new File(inDirS).listFiles();
        List<File> fList = new ArrayList<>();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            if(fs[i].isDirectory()) continue;
            fList.add(fs[i]);
        }
        Collections.sort(fList);
        File[] fsArray = fList.toArray(new File[fList.size()]);
        Arrays.sort(fsArray);
        return fsArray;
    }

    /**
     * 获取文件集合
     * @param inDirS
     * @return
     */
    public static List<File> getFileEndwithInDir (String inDirS, String suffix) {
        File[] fs = new File(inDirS).listFiles();
        List<File> fList = new ArrayList<>();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            if(fs[i].isDirectory()) continue;
            if (!fs[i].getName().endsWith(suffix)) continue;
            fList.add(fs[i]);
        }
        Collections.sort(fList);
        return fList;
    }

    /**
     * 获取文件集合
     * @param inDirS
     * @return
     */
    public static List<File> getFileListInDir (String inDirS) {
        File[] fs = new File(inDirS).listFiles();
        List<File> fList = new ArrayList<>();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            if(fs[i].isDirectory()) continue;
            fList.add(fs[i]);
        }
        Collections.sort(fList);
        return fList;
    }


    /**
     * ############## 该程序待验证
     * @param fs
     * @param colomnIndexChr
     * @param columnIndexPos
     * @param outfileS
     */
    public static void mergeTxtandChangeChrPos(File[] fs, int colomnIndexChr, int columnIndexPos,String outfileS){
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/001_depth/004_50000_Sites/004_AddReliableIntersectGroup";
//        File[] fs = new File(infileDirS).listFiles();

//        fs = IOUtils.listFilesEndsWith(fs,"_AB_sample.txt.gz");
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/001_depth/004_50000_Sites/005_changeChrPos/AB_Popdepth_sample.txt.gz";

//        fs = IOUtils.listFilesEndsWith(fs,"_ABD_sample.txt.gz");
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/001_depth/004_50000_Sites/005_changeChrPos/ABD_Popdepth_sample.txt.gz";

//        fs = IOUtils.listFilesEndsWith(fs,"_D_sample.txt.gz");
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/001_depth/004_50000_Sites/005_changeChrPos/D_Popdepth_sample.txt.gz";

        try{
            Arrays.sort(fs);
            String infileS = fs[0].getAbsolutePath();
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            //read header
            bw.write( br.readLine());
            bw.newLine();

            int cnttotal=0;
            for (int i = 0; i < fs.length; i++) {
                infileS = fs[i].getAbsolutePath();

                br = AoFile.readFile(infileS);
                br.readLine();
                String temp = null; //read header
                int cnt = 0;
                List<String> l = new ArrayList<>();
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    cnt++;
                    cnttotal++;
                    int chrID = Integer.parseInt(l.get(colomnIndexChr));
                    int pos = Integer.parseInt(l.get(columnIndexPos));
                    String chr = RefV1Utils.getChromosome(chrID,pos);
                    int posOnChrosome = RefV1Utils.getPosOnChromosome(chrID,pos);
                    //先找到 chr 所在的列
                    StringBuilder sb = new StringBuilder();
                    for (int j = 0; j < l.size(); j++) {
                        if (j == colomnIndexChr){
                            sb.append(chr).append("\t");
                        }
                        if (j == columnIndexPos){
                            sb.append(posOnChrosome).append("\t");
                        }
                        else{
                            sb.append(l.get(j)).append("\t");
                        }
                    }
                    sb.deleteCharAt(sb.length()-1);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                System.out.println(fs[i].getName() + "\t" + cnt);
            }
            System.out.println("Total lines without header count is " + cnttotal + " at merged file " + outfileS );
            br.close();
            bw.flush();
            bw.close();

        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    /**
     *
     * @param outfileS
     */
    public static void mergeTxtbysuffix(File[] fs, String outfileS, String suffix) {
        fs = IOUtils.listFilesContains(fs, suffix);
        Arrays.sort(fs);
        try {
            String infileS = fs[0].getAbsolutePath();
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            //read header
            bw.write(br.readLine());
            bw.newLine();

            int cnttotal = 0;
            //read context
            for (int i = 0; i < fs.length; i++) {
                infileS = fs[i].getAbsolutePath();
                br = AoFile.readFile(infileS);
                br.readLine();
                String temp = null; //read header
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    cnttotal++;
                    StringBuilder sb = new StringBuilder();
                    sb.append(temp);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                System.out.println(fs[i].getName() + "\t" + cnt);
            }
            System.out.println("Total lines without header count is " + cnttotal + " at merged file " + outfileS );
            br.close();
            bw.flush();
            bw.close();

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     *
     * @param infileDirS
     * @param outfileS
     */
    public static void mergeTxtbysuffix(String infileDirS, String outfileS, String suffix) {
        File[] fs = AoFile.getFileArrayInDir(infileDirS);
        fs = IOUtils.listFilesContains(fs, suffix);
        Arrays.sort(fs);
        try {
            String infileS = fs[0].getAbsolutePath();
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            //read header
            bw.write(br.readLine());
            bw.newLine();

            int cnttotal = 0;
            //read context
            for (int i = 0; i < fs.length; i++) {
                infileS = fs[i].getAbsolutePath();
                br = AoFile.readFile(infileS);
                br.readLine();
                String temp = null; //read header
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    cnttotal++;
                    StringBuilder sb = new StringBuilder();
                    sb.append(temp);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                System.out.println(fs[i].getName() + "\t" + cnt);
            }
            System.out.println("Total lines without header count is " + cnttotal + " at merged file " + outfileS );
            br.close();
            bw.flush();
            bw.close();

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }


    public static void mergeTxt_byFileArray(File[] fs, String outfileS) {
        Arrays.sort(fs);
        try {
            String infileS = fs[0].getAbsolutePath();
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            //read header
            bw.write(br.readLine());
            bw.newLine();

            int cnttotal = 0;
            //read context
            for (int i = 0; i < fs.length; i++) {
                infileS = fs[i].getAbsolutePath();
                br = AoFile.readFile(infileS);
                br.readLine(); //read header
                String temp = null;
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    cnttotal++;
                    StringBuilder sb = new StringBuilder();
                    sb.append(temp);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                System.out.println(fs[i].getName() + "\t" + cnt);
            }
            System.out.println("Total lines without header count is " + cnttotal + " at merged file " + outfileS );
            br.close();
            bw.flush();
            bw.close();

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }



    /**
     * 合并m没有表头的文件,并手动添加表头
     * @param infileDirS
     * @param outfileS
     */
    public static void mergeTxtwithoutHeader(String infileDirS, String header, String outfileS) {
        File[] fs = AoFile.getFileArrayInDir(infileDirS);
        Arrays.sort(fs);
        try {
            BufferedReader br = null;
            BufferedWriter bw = AoFile.writeFile(outfileS);
            //read header
            bw.write(header);
            bw.newLine();
            int cnttotal = 0;
            //read context
            for (int i = 0; i < fs.length; i++) {
                String infileS = fs[i].getAbsolutePath();
                br = AoFile.readFile(infileS);
//                br.readLine();
                String temp = null; //read header
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    cnttotal++;
                    StringBuilder sb = new StringBuilder();
                    sb.append(temp);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                System.out.println(fs[i].getName() + "\t" + cnt);
            }
            System.out.println("Total lines without header count is " + cnttotal + " at merged file " + outfileS );
            br.close();
            bw.flush();
            bw.close();

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }


    /**
     *
     * @param infileDirS
     * @param outfileS
     */
    public static void mergeTxt(String infileDirS, String outfileS) {
        File[] fs = AoFile.getFileArrayInDir(infileDirS);
        Arrays.sort(fs);
        try {
            String infileS = fs[0].getAbsolutePath();
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            //read header
            bw.write(br.readLine());
            bw.newLine();

            int cnttotal = 0;
            //read context
            for (int i = 0; i < fs.length; i++) {
                infileS = fs[i].getAbsolutePath();
                br = AoFile.readFile(infileS);
                br.readLine();
                String temp = null; //read header
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    cnttotal++;
                    StringBuilder sb = new StringBuilder();
                    sb.append(temp);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                System.out.println(fs[i].getName() + "\t" + cnt);
            }
            System.out.println("Total lines without header count is " + cnttotal + " at merged file " + outfileS );
            br.close();
            bw.flush();
            bw.close();

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 目的：将所有txt文本的chr pos位点等信息合并成一个文件。
     *
     * @param infileDirS
     * @param outfileS
     */
    public static void mergeTxtwithoutHeader(String infileDirS, String outfileS) {

        File[] fs = AoFile.getFileArrayInDir(infileDirS);
        Arrays.sort(fs);
        try {
            BufferedReader br = null;
            BufferedWriter bw = AoFile.writeFile(outfileS);
            int cnttotal = 0;
            //读正文部分
            for (int i = 0; i < fs.length; i++) {
                String infileS = fs[i].getAbsolutePath();
                br = AoFile.readFile(infileS);
                String temp = null;
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    cnttotal++;
                    StringBuilder sb = new StringBuilder();
                    sb.append(temp);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                System.out.println(fs[i].getName() + "\t" + cnt);
            }
            System.out.println("Total lines without header count is " + cnttotal + " at merged file " + outfileS );
            br.close();
            bw.flush();
            bw.close();

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }


    public static BufferedReader readFile(String infileS){
        BufferedReader br = null;
        if (infileS.endsWith(".vcf")) {
            br = IOUtils.getTextReader(infileS);
        } else if (infileS.endsWith(".vcf.gz")) {
            br = IOUtils.getTextGzipReader(infileS);
        }
        if (infileS.endsWith(".txt")) {
            br = IOUtils.getTextReader(infileS);
        } else if (infileS.endsWith(".txt.gz")) {
            br = IOUtils.getTextGzipReader(infileS);
        }else if (infileS.endsWith(".hwe.gz")) {
            br = IOUtils.getTextGzipReader(infileS);
        }else if (infileS.endsWith(".csv.gz")) {
            br = IOUtils.getTextGzipReader(infileS);
        }else if (infileS.endsWith(".gz")) {
            br = IOUtils.getTextGzipReader(infileS);
        }else if (infileS.endsWith(".TajimaD.gz")) {
            br = IOUtils.getTextGzipReader(infileS);
        }
        else if (infileS.endsWith(".fst")) {
            br = IOUtils.getTextReader(infileS);
        }else if (infileS.endsWith(".Tajima.D")) {
            br = IOUtils.getTextReader(infileS);
        }
        else if (infileS.endsWith(".pi")) {
            br = IOUtils.getTextReader(infileS);
        }
        else if (infileS.endsWith(".TAB")) {
            br = IOUtils.getTextReader(infileS);
        }
        else if (infileS.endsWith(".csv")) {
            br = IOUtils.getTextReader(infileS);
        }
        else if (infileS.endsWith(".tsv")) {
            br = IOUtils.getTextReader(infileS);
        }
        else if (infileS.endsWith(".sh")) {
            br = IOUtils.getTextReader(infileS);
        }
        else if (infileS.endsWith(".md5")) {
            br = IOUtils.getTextReader(infileS);
        }
        else if (infileS.endsWith(".hwe")) {
            br = IOUtils.getTextReader(infileS);
        }
        else if (infileS.endsWith(".fa")) {
            br = IOUtils.getTextReader(infileS);
        }else if (infileS.endsWith(".fasta")) {
            br = IOUtils.getTextReader(infileS);
        }
        return br;
    }


    public static void writeListtoFile(List<String> list, String outfileS){
        try {
            BufferedWriter bw = AoFile.writeFile(outfileS);
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < list.size(); i++) {
                sb.setLength(0);
                sb.append(list.get(i));
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public static BufferedWriter writeFile(String outfileS){
        BufferedWriter bw = null;
        if (outfileS.endsWith(".txt")){
            bw=IOUtils.getTextWriter(outfileS);
        }else if (outfileS.endsWith(".txt.gz")) {
            bw = IOUtils.getTextGzipWriter(outfileS);
        }
        if (outfileS.endsWith(".vcf")){
            bw=IOUtils.getTextWriter(outfileS);
        }else if (outfileS.endsWith(".vcf.gz")) {
            bw = IOUtils.getTextGzipWriter(outfileS);
        }else if (outfileS.endsWith(".csv.gz")) {
            bw = IOUtils.getTextGzipWriter(outfileS);
        }else if (outfileS.endsWith(".sh")){
            bw=IOUtils.getTextWriter(outfileS);
        }else if (outfileS.endsWith(".md5")){
            bw=IOUtils.getTextWriter(outfileS);
        }else if (outfileS.endsWith(".fst")){
            bw=IOUtils.getTextWriter(outfileS);
        }else if (outfileS.endsWith(".hwe")) {
            bw = IOUtils.getTextWriter(outfileS);
        }else if (outfileS.endsWith(".csv")) {
            bw = IOUtils.getTextWriter(outfileS);
        }else if (outfileS.endsWith(".fasta")) {
            bw = IOUtils.getTextWriter(outfileS);
        }
        return bw;
    }

    public static void addGroupbyFileName(String infileS, String regularE, String outfileS){


    }


    /**
     * add more than one column to a file
     */
    public static void addColumsbyString(String infileS,int keyIDindex, HashMap<String,String>[] hm,String headername){
        String outfileS = null;
        String outfileDirS = new File(infileS).getParent(); //获取输入文件的目录
        outfileDirS = new File(outfileDirS).getParent(); //根据输入文件的目录获取上一级父目录；
        outfileDirS = outfileDirS + "/A_out"; //根据父目录路径 创建输出文件目录
        new File(outfileDirS).mkdirs(); //创建输出文件目录
        outfileS = new File(outfileDirS,new File(infileS).getName()).getAbsolutePath();

        try{
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);

            String temp = br.readLine(); //read header
            bw.write(temp + "\t" + headername);
            bw.newLine();
            List<String> l = new ArrayList<>();
            while((temp=br.readLine()) != null){
                l = PStringUtils.fastSplit(temp);
                String key = l.get(keyIDindex); //注意，如果string类型不能转化为pos,这里也不会报错
                bw.write(temp);

                for (int i = 0; i < hm.length; i++) {
                    String value = hm[i].get(key);
                    if(value == null || value == ""){  //!!!!! if there is no value, we should set the value as "NA".
                        value = "NA";
                    }
                    bw.write("\t" + value);
                }
                bw.newLine();
            }
            bw.flush();
            bw.close();
            br.close();
            System.out.println(infileS + "\tis completed at\t" + outfileS);
        }
        catch(Exception e){
            System.exit(1);
            e.printStackTrace();
        }
    }


    /**
     * add colum to a file
     */
    public static void addColumbyString(String infileS,int keyIDindex, HashMap<String,String> hm,String headername){
        String outfileS = null;
        String outfileDirS = new File(infileS).getParent(); //获取输入文件的目录
        outfileDirS = new File(outfileDirS).getParent(); //根据输入文件的目录获取上一级父目录；
        outfileDirS = outfileDirS + "/A_out"; //根据父目录路径 创建输出文件目录
        new File(outfileDirS).mkdirs(); //创建输出文件目录
        outfileS = new File(outfileDirS,new File(infileS).getName()).getAbsolutePath();

        try{
            BufferedReader br = null;
            BufferedWriter bw = null;
            if(infileS.endsWith(".txt")){
                br=IOUtils.getTextReader(infileS);
                bw = IOUtils.getTextWriter(outfileS);
            }else if(infileS.endsWith(".txt.gz")){
                br=IOUtils.getTextGzipReader(infileS);
                bw=IOUtils.getTextGzipWriter(outfileS);
            }

            String temp = br.readLine(); //read header
            bw.write(temp + "\t" + headername);
            bw.newLine();
            List<String> l = new ArrayList<>();
            while((temp=br.readLine()) != null){
                l = PStringUtils.fastSplit(temp);
                String key = l.get(keyIDindex); //注意，如果string类型不能转化为pos,这里也不会报错
                String value = hm.get(key);
                if(value == null || value == ""){  //!!!!! if there is no value, we should set the value as "NA".
                    value = "NA";
                }
                bw.write(temp + "\t" + value);
                bw.newLine();
            }
            bw.flush();
            bw.close();
            br.close();
            System.out.println(infileS + "\tis completed at\t" + outfileS);
        }
        catch(Exception e){
            System.exit(1);
            e.printStackTrace();
        }
    }


    /**
     * add colum to a file
     * 注意： header的第一个单词不用加制表符
     */
    public static void addColumbyint(String infileS,int keyIDindex, HashMap<Integer,String> hm,String headername){
        String outfileS = null;
        String outfileDirS = new File(infileS).getParent(); //获取输入文件的父目录
        outfileDirS = new File(outfileDirS).getParent(); //根据输入文件的父目录获取上一级父目录；
        outfileDirS = outfileDirS + "/A_out"; //根据父目录路径 创建输出文件目录
        new File(outfileDirS).mkdirs(); //创建输出文件目录
        outfileS = new File(outfileDirS,new File(infileS).getName()).getAbsolutePath();

        try{
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = br.readLine(); //read header
            bw.write(temp + "\t" + headername);
            bw.newLine();
            List<String> l = new ArrayList<>();
            while((temp=br.readLine()) != null){
                l = PStringUtils.fastSplit(temp);
                int key = Integer.parseInt(l.get(keyIDindex)); //注意，如果string类型不能转化为pos,这里也不会报错
                String value = hm.get(key);
                if(value == null || value == "" || value.isEmpty()){  //!!!!! if there is no value, we should set the value as "NA".
                    value = "NA";
                }
                bw.write(temp + "\t" + value);
                bw.newLine();
            }
            bw.flush();
            bw.close();
            br.close();
            System.out.println(infileS + "\tis completed at\t" + outfileS);
        }
        catch(Exception e){
            System.exit(1);
            e.printStackTrace();
        }
    }


    /**
     * return a hashmap from a file, the value is the file name
     *
     * @param infileS
     * @return
     */
    public HashMap<String,String> getHashMapwithFilename(String infileS){
        HashMap<String,String> hm = new HashMap<>();
        try{
            String group = new File(infileS).getName().replaceFirst(".txt","");
            BufferedReader br = IOUtils.getTextReader(infileS);
            String temp = null; //do not read header
            int cnt = 0;
            while((temp = br.readLine()) != null){
                cnt++;
                hm.put(temp,group);
            }
            br.close();
            System.out.println(new File(infileS).getName() + "\tHashMap size is\t" + cnt);
            System.out.println("Total HashMap size is " + cnt);
        }
        catch(Exception e){
            System.exit(1);
            e.printStackTrace();
        }

        System.out.println("HashMap contains " + hm.size() + " pairs");

        return hm;
    }

    /**
     * return a hashmap from a file dirs
     * 每个文件的名字是分组信息，文件内的taxa是key
     * 此目录下的所有文件都建立这样的关系，返回到同一个hashmap中
     *
     * @param infileDirS
     * @return
     */
    public HashMap<String,String> getHashMapwithFileDirs(String infileDirS){
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                System.out.println(fs[i].getName() + " is hidden");
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        Arrays.sort(fs);
        int cnttotal = 0;
        HashMap<String,String> hm = new HashMap<>();
        try{
            for (int i = 0; i < fs.length; i++) {
                String infileS = fs[i].getAbsolutePath();
                String group = fs[i].getName().replaceFirst(".txt","");
                BufferedReader br = IOUtils.getTextReader(infileS);
                String temp = null; //do not read header
                int cnt = 0;
                while((temp = br.readLine()) != null){
                    cnttotal++;
                    cnt++;
                    hm.put(temp,group);
                }
                br.close();
                System.out.println(fs[i].getName() + "\tHashMap size is\t" + cnt);
            }
            System.out.println("Total HashMap size is " + cnttotal);
        }
        catch(Exception e){
            System.exit(1);
            e.printStackTrace();
        }

        System.out.println("HashMap contains " + hm.size() + " pairs");


        return hm;
    }


    /**
     *return a hashmap from a file
     *
     * @param infileS
     * @param keycolummIndex
     * @param valuecolumnIndex
     * @return
     */
    public static HashMap<String,Double> getHashMapdoubleValue(String infileS, int keycolummIndex, int valuecolumnIndex){


        String out = null;
        RowTable<String> t = new RowTable<>(infileS);
        HashMap<String,Double> hm = new HashMap<>();
        for (int i = 0; i < t.getRowNumber() ; i++) {
            String key = t.getCell(i,keycolummIndex);
            String value = t.getCell(i,valuecolumnIndex);
            double valued = Double.parseDouble(value);
            hm.put(key,valued);
        }
        System.out.println("HashMap contains " + hm.size() + " pairs");


        return hm;
    }

    /**
     *return a hashmap from a file
     *
     * @param infileS
     * @param keycolummIndex
     * @param valuecolumnIndex
     * @return
     */
    public static HashMap<Integer,String> getHashMapintKey(String infileS, int keycolummIndex, int valuecolumnIndex){

        String out = null;
        RowTable<String> t = new RowTable<>(infileS);
        HashMap<Integer,String> hm = new HashMap<>();
        for (int i = 0; i < t.getRowNumber() ; i++) {
            int key = Integer.parseInt(t.getCell(i,keycolummIndex));
            String value = t.getCell(i,valuecolumnIndex);
            hm.put(key,value);
        }
        System.out.println("HashMap contains " + hm.size() + " pairs");


        return hm;
    }

    /**
     *return a hashmap from a file
     *
     * @param infileS
     * @param keycolummIndex
     * @param valuecolumnIndexs
     * @return
     */
    public static HashMap<String,String>[] getHashMapsStringKey(String infileS, int keycolummIndex, int[] valuecolumnIndexs){
        String out = null;
        RowTable<String> t = new RowTable<>(infileS);
        HashMap<String,String>[] hm = new HashMap[valuecolumnIndexs.length];
        for (int i = 0; i < hm.length; i++) {
            hm[i] = new HashMap<String,String>();
        }
        for (int i = 0; i < t.getRowNumber() ; i++) {
            String key = t.getCell(i,keycolummIndex);
            for (int j = 0; j < valuecolumnIndexs.length; j++) {
                String value = t.getCell(i,valuecolumnIndexs[j]);
                hm[j].put(key,value);
            }
        }
        System.out.println("HashMap contains " + hm[0].size() + " pairs and totally " + valuecolumnIndexs.length + " HashMaps are built");


        return hm;
    }

    /**
     *return a hashmap from a file
     *
     * @param infileS
     * @param keycolummIndex
     * @param valuecolumnIndex
     * @return
     */
    public static HashMap<String,String> getHashMapStringKey(String infileS, int keycolummIndex, int valuecolumnIndex){
        String out = null;
        RowTable<String> t = new RowTable<>(infileS);
        HashMap<String,String> hm = new HashMap<>();
        for (int i = 0; i < t.getRowNumber() ; i++) {
            String key = t.getCell(i,keycolummIndex);
            String value = t.getCell(i,valuecolumnIndex);
            hm.put(key,value);
        }
        System.out.println("HashMap contains " + hm.size() + " pairs");

        
        return hm;
    }

    /**
     *return a hashmap from a file
     *
     * @param infileS
     * @param keycolummIndex
     * @param valuecolumnIndex
     * @return
     */
    public static HashMap<String,String> getHashMapStringKey_withoutHeader(String infileS, int keycolummIndex, int valuecolumnIndex){
        HashMap<String,String> hm = new HashMap<>();
        try {
            BufferedReader br = AoFile.readFile(infileS);
            String temp = null; //no header
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                String key = l.get(keycolummIndex);
                String value = l.get(valuecolumnIndex);
                hm.put(key,value);
            }
            br.close();
            System.out.println("The file has " + cnt + " lines without counting header.");
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        System.out.println("HashMap contains " + hm.size() + " pairs");


        return hm;
    }

//######################################
    //##################################
    public class Record implements Comparable<Record>{
        public String taxa;
        public String sub;

        public Record(String taxa, String sub) {
            this.taxa = taxa;
            this.sub = sub;

        }
        public boolean isSimilar(int pos, String alt) {
            if (taxa.equals(this.taxa) && sub.equals(this.sub)) {
                return true;
            }
            return false;
        }

        @Override
        public int compareTo(Record o) {
            if (this.taxa.compareTo(o.taxa)< 0){
                return -1;

            }else if (this.taxa.compareTo(o.taxa)==0){
                return this.sub.compareTo(o.sub);
            }else{
                return 1;
            }
        }
    }

    /**
     *
     * get HashMap list from a txt file
     * @param infileS
     * @param columnIndex1
     * @param columnIndex2
     * @return
     */
    public List<Record> getRecordList(String infileS, int columnIndex1, int columnIndex2){
        List<Record> out = new ArrayList<>();
        try {
            BufferedReader br = null;
            if (infileS.endsWith(".txt")) {
                br = IOUtils.getTextReader(infileS);
            } else if (infileS.endsWith(".txt.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }

            String temp = br.readLine(); //read header
            List<String> l = new ArrayList();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String goal1 = l.get(columnIndex1);
                String goal2 = l.get(columnIndex2);
                Record r = new Record(goal1,goal2);
                out.add(r);
                cnt++;
            }
            br.close();
            System.out.println("Total num in the list is    " + cnt + "\t" + out.size());
            for (int i = 0; i < out.size(); i++) {
                System.out.println(out.get(i));
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return out;
    }

    //##################################
    //##################################

    /**
     *
     * get HashMap list from a txt file
     * @param infileS
     * @param columnIndex1
     * @param columnIndex2
     * @return
     */
    public List<HashMap<String,String>> getHashMapList(String infileS, int columnIndex1, int columnIndex2){
        List<HashMap<String,String>> out = new ArrayList<>();
        try {
            BufferedReader br = null;
            if (infileS.endsWith(".txt")) {
                br = IOUtils.getTextReader(infileS);
            } else if (infileS.endsWith(".txt.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }

            String temp = br.readLine(); //read header
            List<String> l = new ArrayList();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String goal1 = l.get(columnIndex1);
                String goal2 = l.get(columnIndex2);
                HashMap<String,String> hm = new HashMap<>();
                hm.put(goal1,goal2);
                out.add(hm);
                cnt++;
            }
            br.close();
            System.out.println("Total num in the list is    " + cnt + "\t" + out.size());
            for (int i = 0; i < out.size(); i++) {
                System.out.println(out.get(i));
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return out;
    }

    /**
     *
     * get String list from a txt file which contains group set
     * @param infileS
     * @param columnIndex
     * @return
     */
    public static List<String> getStringListbySet(String infileS, int columnIndex){

        Set<String> out = new HashSet<>();
        try {
            BufferedReader br = null;
            if (infileS.endsWith(".txt")) {
                br = IOUtils.getTextReader(infileS);
            } else if (infileS.endsWith(".txt.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }

            String temp = br.readLine(); //read header
            List<String> l = new ArrayList();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String goal = l.get(columnIndex);
                out.add(goal);
                cnt++;
            }
            br.close();
            System.out.println("Total num in the set is " + out.size());
        }
        catch (Exception e) {
            e.printStackTrace();
        }

        List<String> out2 = new ArrayList<>(out);
        Collections.sort(out2);
        return out2;
    }

    /**
     *
     * get String set from a txt file
     * @param infileS
     * @param columnIndex
     * @return
     */
    public static String[] getStringArraybySet(String infileS, int columnIndex){
        Set<String> out = new HashSet<>();
        try {
            BufferedReader br = AoFile.readFile(infileS);
            String temp = br.readLine(); //read header
            List<String> l = new ArrayList();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String goal = l.get(columnIndex);
                out.add(goal);
                cnt++;
            }
            br.close();
            System.out.println("Total num in the set is " + out.size());
        }
        catch (Exception e) {
            e.printStackTrace();
        }

        String[] outArray = out.toArray(new String[out.size()]);
        Arrays.sort(outArray);
        return outArray;
    }

    /**
     *
     * get String set from a txt file
     * @param infileS
     * @param columnIndex
     * @return
     */
    public Set<String> getStringSet(String infileS, int columnIndex){
        Set<String> out = new HashSet<>();
        try {
            BufferedReader br = null;
            if (infileS.endsWith(".txt")) {
                br = IOUtils.getTextReader(infileS);
            } else if (infileS.endsWith(".txt.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }

            String temp = br.readLine(); //read header
            List<String> l = new ArrayList();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String goal = l.get(columnIndex);
                out.add(goal);
                cnt++;
            }
            br.close();
            System.out.println("Total num in the set is " + out.size());
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return out;
    }


    /**
     * get String list from a txt file
     *
     * @param infileS
     * @param columnIndex
     * @return
     */
    public static List<String> getStringListwithoutHeader(String infileS, int columnIndex){
        List<String> out = new ArrayList<>();
        try {
            BufferedReader br = null;
            if (infileS.endsWith(".txt")) {
                br = IOUtils.getTextReader(infileS);
            } else if (infileS.endsWith(".txt.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }

            String temp = null; //read header
            List<String> l = new ArrayList();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String goal = l.get(columnIndex);
                out.add(goal);
                cnt++;
            }
            br.close();
            System.out.println("Total num in the list is    " + cnt + "\t" + out.size());
            Collections.sort(out);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return out;
    }

    /**
     *
     * get String list from a txt file
     * @param infileS
     * @param columnIndex
     * @return
     */
    public static List<String> getStringList(String infileS, int columnIndex){
        List<String> out = new ArrayList<>();
        try {
            BufferedReader br = null;
            if (infileS.endsWith(".txt")) {
                br = IOUtils.getTextReader(infileS);
            } else if (infileS.endsWith(".txt.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }

            String temp = br.readLine(); //read header
            List<String> l = new ArrayList();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String goal = l.get(columnIndex);
                out.add(goal);
                cnt++;
            }
            br.close();
            System.out.println("Total num in the list is"+ "\t" + out.size());
            Collections.sort(out);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return out;
    }

    /**
     *
     * get String list from a txt file
     * @param infileS
     * @param columnIndex
     * @return
     */
    public static String[] getStringArraybyList_withoutHeader(String infileS, int columnIndex){

        List<String> out = new ArrayList<>();
        try {
            BufferedReader br = null;
            if (infileS.endsWith(".txt")) {
                br = IOUtils.getTextReader(infileS);
            } else if (infileS.endsWith(".txt.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }

            String temp = null; //read header
            List<String> l = new ArrayList();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String goal = l.get(columnIndex);
                out.add(goal);
                cnt++;
            }
            br.close();
            System.out.println("Total num in the list is    " + cnt + "\t" + out.size());
            Collections.sort(out);
        }
        catch (Exception e) {
            e.printStackTrace();
        }

        String[] outArray = out.toArray(new String[out.size()]);
        Arrays.sort(outArray);
        return outArray;
    }

    /**
     *
     * get String list from a txt file
     * @param infileS
     * @param columnIndex
     * @return
     */
    public static String[] getgeneArraybyList(String infileS, int columnIndex){

        List<String> out = new ArrayList<>();
        try {
            BufferedReader br = AoFile.readFile(infileS);
            String temp = br.readLine(); //read header
            List<String> l = new ArrayList();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String goal = l.get(columnIndex);
                out.add(goal.split("\\.")[0]);
                cnt++;
            }
            br.close();
            System.out.println("Total num in the list is    " + cnt + "\t" + out.size());
            Collections.sort(out);
        }
        catch (Exception e) {
            e.printStackTrace();
        }

        String[] outArray = out.toArray(new String[out.size()]);
        Arrays.sort(outArray);
        return outArray;
    }

    /**
     *
     * get String list from a txt file
     * @param infileS
     * @param columnIndex
     * @return
     */
    public static String[] getStringArraybyList(String infileS, int columnIndex){

        List<String> out = new ArrayList<>();
        try {
            BufferedReader br = null;
            if (infileS.endsWith(".txt")) {
                br = IOUtils.getTextReader(infileS);
            } else if (infileS.endsWith(".txt.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }

            String temp = br.readLine(); //read header
            List<String> l = new ArrayList();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String goal = l.get(columnIndex);
                out.add(goal);
                cnt++;
            }
            br.close();
            System.out.println("Total num in the list is    " + cnt + "\t" + out.size());
            Collections.sort(out);
        }
        catch (Exception e) {
            e.printStackTrace();
        }

        String[] outArray = out.toArray(new String[out.size()]);
        Arrays.sort(outArray);
        return outArray;
    }


    /**
     * get the pos database from txt file without header
     * @param infileS
     * @param columnIndex
     * @return
     */
    public static TIntArrayList getTIntList_withoutHeader(String infileS, int columnIndex){
        TIntArrayList ll = new TIntArrayList();
        try {
            BufferedReader br = AoFile.readFile(infileS);
            String temp = null;
            List<String> l = new ArrayList();
            int cnttotal = 0;
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                cnttotal++;
                l = PStringUtils.fastSplit(temp);
                String goal = l.get(columnIndex);
                if (goal.startsWith("N") | goal.startsWith("inf")) continue;
                ll.add(Integer.parseInt(goal));
                cnt++;
            }
            br.close();
            System.out.println("Total num in the list is " + cnt + "\tTDoubleArrayList size is " + ll.size());
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return ll;
    }

    /**
     * get the pos database from txt file
     *
     * @param infileS
     * @param columnIndex
     * @return
     */
    public static TIntArrayList getTIntList(String infileS, int columnIndex){
        TIntArrayList ll = new TIntArrayList();
        try {
            BufferedReader br = AoFile.readFile(infileS);
            String temp = br.readLine(); //read header
            List<String> l = new ArrayList();
            int cnttotal = 0;
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                cnttotal++;
                l = PStringUtils.fastSplit(temp);
                String goal = l.get(columnIndex);
                if (goal.startsWith("N") | goal.startsWith("inf")) continue;
                ll.add(Integer.parseInt(goal));
                cnt++;
            }
            br.close();
            System.out.println("Total num in the list is " + cnt + "\tTDoubleArrayList size is " + ll.size());
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return ll;
    }

    /**
     * get the pos database from txt file
     *
     * @param infileS
     * @param columnIndex
     * @return
     */
    public static TDoubleArrayList getTDoubleList(String infileS, int columnIndex){
        TDoubleArrayList ll = new TDoubleArrayList();
        try {
            BufferedReader br = AoFile.readFile(infileS);
            String temp = br.readLine(); //read header
            List<String> l = new ArrayList();
            int cnttotal = 0;
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                cnttotal++;
                l = PStringUtils.fastSplit(temp);
                String goal = l.get(columnIndex);
                if(goal == null || goal == "" || goal.isEmpty()){  //!!!!! if there is no value, we should set the value as "NA".
                    goal = "NA";
                }
                if (goal.startsWith("N") || goal.startsWith("inf") || goal.startsWith("-") || goal.startsWith("I")) continue; //
                ll.add(Double.parseDouble(goal));
                cnt++;
            }
            br.close();
            System.out.println("Total num in the list is " + cnt + "\tTDoubleArrayList size is " + ll.size());
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return ll;
    }

    /**
     * get the pos database from txt file
     *
     * @param infileS
     * @param columnIndex
     * @return
     */
    public TIntArrayList getNumList(String infileS, int columnIndex){
        TIntArrayList ll = new TIntArrayList();

        try {
            BufferedReader br = AoFile.readFile(infileS);
            String temp = br.readLine(); //read header
            List<String> l = new ArrayList();
            int cnttotal = 0;
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                    cnttotal++;
                    l = PStringUtils.fastSplit(temp);
                    String goal = l.get(columnIndex);
                    if (goal.startsWith("N")) continue;
                    ll.add(Integer.parseInt(goal));
                    cnt++;
            }
            br.close();
            System.out.println("Total num in the list is    " + cnt + "\t" + ll.size());
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return ll;
    }


    /**
     * get the pos database from VCF file
     *
     * @param infileS
     * @return
     */
    public static TIntArrayList getNumListfromVCF(String infileS){
        TIntArrayList ll = new TIntArrayList();

        try {
            BufferedReader br = AoFile.readFile(infileS);
            String temp = null;
            List<String> l = new ArrayList();
            int cnttotal = 0;
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                if (!temp.startsWith("#")) {
                    cnttotal++;
                    l = PStringUtils.fastSplit(temp);
                    String goal = l.get(1);
                    if (goal.startsWith("N")) continue;
                    ll.add(Integer.parseInt(goal));
                    cnt++;
                }
            }
            br.close();
            System.out.println("Total num in the list is    " + cnt + "\t" + ll.size());
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return ll;
    }

    /**
     *
     * @param infileS
     */
    public static List<String> getheader(String infileS){

        List<String> l = new ArrayList<>();
        try {
            BufferedReader br = null;
            if (infileS.endsWith(".txt")) {
                br = IOUtils.getTextReader(infileS);
            } else if (infileS.endsWith(".txt.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }
            else if (infileS.endsWith(".vcf.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }else if (infileS.endsWith(".vcf")) {
                br = IOUtils.getTextReader(infileS);
            }
            else if (infileS.endsWith(".xls")) {
                br = IOUtils.getTextReader(infileS);
            }
            else if (infileS.endsWith(".TAB")) {
                br = IOUtils.getTextReader(infileS);
            }
            else if (infileS.endsWith(".tsv")) {
                br = IOUtils.getTextReader(infileS);
            }
            else if (infileS.endsWith(".gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }
            String temp = br.readLine();
            l = PStringUtils.fastSplit(temp);
            br.close();

        } catch (Exception e) {
            e.printStackTrace();
        }

        return l;

    }

    /**
     *
     * @param infileS
     */
    public static List<String> getHeader(String infileS){
        List<String> l = new ArrayList<>();
        try {
            BufferedReader br = null;
            if (infileS.endsWith(".txt")) {
                br = IOUtils.getTextReader(infileS);
            } else if (infileS.endsWith(".txt.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }
            else if (infileS.endsWith(".vcf.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }else if (infileS.endsWith(".vcf")) {
                br = IOUtils.getTextReader(infileS);
            }
            else if (infileS.endsWith(".xls")) {
                br = IOUtils.getTextReader(infileS);
            }
            else if (infileS.endsWith(".TAB")) {
                br = IOUtils.getTextReader(infileS);
            }
            else if (infileS.endsWith(".tsv")) {
                br = IOUtils.getTextReader(infileS);
            }
            else if (infileS.endsWith(".gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }
            String temp = br.readLine();
            l = PStringUtils.fastSplit(temp);
            int cnt = -1;
//            for (int i = 0; i < l.size(); i++) {
//                System.out.println(String.valueOf(i) + "\t" + String.valueOf(l.get(i)));
//            }

            br.close();


        } catch (Exception e) {
            e.printStackTrace();
        }
        return l;

    }


    /**
     *
     * @param infileS
     */
    public static void readheader(String infileS){
        try {
            BufferedReader br = null;
            if (infileS.endsWith(".txt")) {
                br = IOUtils.getTextReader(infileS);
            } else if (infileS.endsWith(".txt.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }
            else if (infileS.endsWith(".vcf.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }else if (infileS.endsWith(".hwe.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }else if (infileS.endsWith(".vcf")) {
                br = IOUtils.getTextReader(infileS);
            }
            else if (infileS.endsWith(".xls")) {
                br = IOUtils.getTextReader(infileS);
            }
            else if (infileS.endsWith(".TAB")) {
                br = IOUtils.getTextReader(infileS);
            }
            else if (infileS.endsWith(".tsv")) {
                br = IOUtils.getTextReader(infileS);
            }
            else if (infileS.endsWith(".gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }else if (infileS.endsWith(".fst")) {
                br = IOUtils.getTextReader(infileS);
            }else if (infileS.endsWith(".Tajima.D")) {
                br = IOUtils.getTextReader(infileS);
            }else if (infileS.endsWith(".pi")) {
                br = IOUtils.getTextReader(infileS);
            } else if (infileS.endsWith(".hwe")) {
                br = IOUtils.getTextReader(infileS);
            }


            String temp = br.readLine();
            List<String> l = PStringUtils.fastSplit(temp);
            int cnt = -1;
            for (int i = 0; i < l.size(); i++) {
                System.out.println(String.valueOf(i) + "\t" + String.valueOf(l.get(i)));
            }

            br.close();

        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    /**
     *
     * @param infileS
     * @param outfileS
     */
    public static void subsetTxt_withoutHeaer_singleStream(String infileS, String outfileS, String ratio) {

        try {

            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);

            String temp = null;
            int cnt = 0;
            System.out.println(new SimpleDateFormat().format(new Date()) + "    program execution.\n");
            long startTime = System.nanoTime();
            Double Ratio = Double.parseDouble(ratio);
            while ((temp = br.readLine()) != null) {
                cnt++;
                double r = Math.random();
                if (r > Ratio) continue;
                //返回带正号的 double 值，该值大于等于 0.0 且小于 1.0。返回值是一个伪随机选择的数，在该范围内（近似）均匀分布
                bw.write(temp);
                bw.newLine();
            }
            bw.flush();
            bw.close();
            br.close();
            long endTime = System.nanoTime();
            float excTime = (float) (endTime - startTime) / 1000000000;
            System.out.println("Execution time: " + String.format("%.2f", excTime) + "s");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }


    /**
     * 复制某一目录下的所有文件夹，不复制文件
     *
     * @param infileDirS
     * @param outfileDirS
     */
    public void copyFileDirS(String infileDirS, String outfileDirS){
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/002_changeChrPos";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/003_mergebySub";
        
        
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            System.out.println(fs[i]);
            new File(outfileDirS, fs[i].getName()).mkdirs();
        }
    }
}
