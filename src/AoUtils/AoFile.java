/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package AoUtils;

import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import pgl.format.table.RowTable;
import pgl.utils.IOUtils;
import pgl.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;

/**
 *
 * @author Aoyue
 */
public class AoFile {
    public AoFile(){
        
    }

    public static int getFileRowNumber(String infileS){
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
    public static List<File> getFileListInDir (String inDirS) {
        File[] fs = new File(inDirS).listFiles();
        List<File> fList = new ArrayList<>();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            fList.add(fs[i]);
        }
        Collections.sort(fList);
        return fList;
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

        File[] fs = new AoFile().getFileArrayInDir(infileDirS);
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
        }else if (infileS.endsWith(".gz")) {
            br = IOUtils.getTextGzipReader(infileS);
        }
        else if (infileS.endsWith(".fst")) {
            br = IOUtils.getTextReader(infileS);
        }else if (infileS.endsWith(".pi")) {
            br = IOUtils.getTextReader(infileS);
        }
        else if (infileS.endsWith(".TAB")) {
            br = IOUtils.getTextReader(infileS);
        }
        return br;
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
        }
        return bw;
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
     */
    public void addColumbyint(String infileS,int keyIDindex, HashMap<Integer,String> hm,String headername){
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
    public HashMap<Integer,String> getHashMapintKey(String infileS, int keycolummIndex, int valuecolumnIndex){
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
     * @param valuecolumnIndex
     * @return
     */
    public HashMap<String,String> getHashMapStringKey(String infileS, int keycolummIndex, int valuecolumnIndex){
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
     * get String set from a txt file
     * @param infileS
     * @param columnIndex
     * @return
     */
    public static String[] getStringArraybySet(String infileS, int columnIndex){
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
    public List<String> getStringList(String infileS, int columnIndex){
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
                if (goal.startsWith("N") | goal.startsWith("inf")) continue;
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
            BufferedReader br = null;
            if (infileS.endsWith(".txt")) {
                br = IOUtils.getTextReader(infileS);
            } else if (infileS.endsWith(".txt.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }
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
     * @param columnIndex
     * @return
     */
    public TIntArrayList getNumListfromVCF(String infileS, int columnIndex){
        TIntArrayList ll = new TIntArrayList();

        try {
            BufferedReader br = null;
            if (infileS.endsWith(".vcf.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }else if(infileS.endsWith(".vcf")) {
                br = IOUtils.getTextReader(infileS);
            }

            String temp = null;
            List<String> l = new ArrayList();
            int cnttotal = 0;
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                if (!temp.startsWith("#")) {
                    cnttotal++;
                    l = PStringUtils.fastSplit(temp);
                    String goal = l.get(columnIndex);
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
            }else if (infileS.endsWith(".vcf")) {
                br = IOUtils.getTextReader(infileS);
            }
            else if (infileS.endsWith(".xls")) {
                br = IOUtils.getTextReader(infileS);
            }
            else if (infileS.endsWith(".TAB")) {
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
