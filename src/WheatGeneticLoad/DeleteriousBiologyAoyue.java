package WheatGeneticLoad;

import AoUtils.AoFile;
import AoUtils.CountSites;
import gnu.trove.list.array.TCharArrayList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.table.RowTable;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;

public class DeleteriousBiologyAoyue {

    public DeleteriousBiologyAoyue() {
        ///Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/001_snp/README.txt you can read and make you clear

//        this.mergeDelgenicSNPAnnotation(); //step1
//        this.countDeleteriousVMapII(); //这是算出所有的总和
//        this.countDeleteriousVMapII_byChr(); //step2

        this.mkDepthOfVMapII();
//        this.mkDepthSummary();
//        this.mergeDepthSummary();

//        this.countDeleteriousVMapIIHighDepth(); //step3
//        this.countDeleteriousVMapIIHighDepth_old();
//        this.mergeFinalfilebySub(); //step4


//        this.mkIndexforBurden();
        this.delRatioVSsynRatio();



    }


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
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return out;
    }


    /**
     * 本方法目的是为了查看有害突变的burden和同义突变的burden相对的比率，为了标准化。
     */

    public void delRatioVSsynRatio(){
        //no change this file: means synonymous deriverd allele burden
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/003_delCount_highDepth/003_/additiveDeleterious_vmap2_highDepth_bysub.txt";

        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/005_VMap2.1DelCount/002_delCount_highDepth/003_/additiveDeleterious_vmap2_highDepth_bysub.txt";
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/005_VMap2.1DelCount/002_delCount_highDepth/006_synonymousGerp1/additiveDeleterious_vmap2_highDepth_bysub.txt";

        //can change
//        String delfileS = "";
//        String outfileS = "";

//        String delfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/003_delCount_highDepth/000_/additiveDeleterious_vmap2_highDepth_bysub.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/004_ratiotest/delVSsynratio.txt";

//        String delfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/003_delCount_highDepth/001_/additiveDeleterious_vmap2_highDepth_bysub.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/004_ratiotest/delDerivedSIFTVSsynratio.txt";

//        String delfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/003_delCount_highDepth/002_/additiveDeleterious_vmap2_highDepth_bysub.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/004_ratiotest/delNonsynonymousVSsynonymous.txt";


        // est-sfs
//        String delfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/005_VMap2.1DelCount/002_delCount_highDepth/000_/additiveDeleterious_vmap2_highDepth_bysub.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/005_VMap2.1DelCount/003_ratiotest/delVSsynratio.txt"; //正常sift计算出来的结果


//        String delfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/005_VMap2.1DelCount/002_delCount_highDepth/001_/additiveDeleterious_vmap2_highDepth_bysub.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/005_VMap2.1DelCount/003_ratiotest/delRefaltSIFTVSsynratio.txt"; //正常sift计算出来的结果
//
        String delfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/005_VMap2.1DelCount/002_delCount_highDepth/002_/additiveDeleterious_vmap2_highDepth_bysub.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/005_VMap2.1DelCount/003_ratiotest/delNonsynonymousVSsynonymous.txt"; //正常sift计算出来的结果


//        String delfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/005_VMap2.1DelCount/002_delCount_highDepth/005_nonsynonymousGerp1/additiveDeleterious_vmap2_highDepth_bysub.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/005_VMap2.1DelCount/003_ratiotest/nonsynonymousGerp1VSsynratio.txt";

        List<Record> taxasubl = this.getRecordList(infileS,0,1);
        Collections.sort(taxasubl);

        //建立 1 based 的 hashmap,根据行数找到对应的 synonymous 的 ratio
        Map<Integer,String> hmratio = new HashMap<>();
        RowTable<String> t = new RowTable<>(infileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            String ratio = t.getCellAsString(i,7);
            hmratio.put(i+1,ratio);
        }

        String taxaGroupFileS = "/Users/Aoyue/project/wheatVMapII/001_germplasm/GermplasmDB/001_toFeiLu/004_wheatVMapII_germplasmInfoAddSubspecies_20191225.txt";
        String taxaGroup2FileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/000_prepareData/001_input/taxaList.txt";
        t = new RowTable (taxaGroupFileS);
        HashMap<String, String> taxaGroupMap = new HashMap();
        HashMap<String, String> taxaSubMap = new HashMap();
        HashMap<String, String> taxaGroupIDMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            taxaGroupMap.put(t.getCellAsString(i, 4), t.getCellAsString(i, 26));
            taxaSubMap.put(t.getCellAsString(i, 4), t.getCellAsString(i, 27));
        }
        t = new RowTable (taxaGroup2FileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            taxaGroupIDMap.put(t.getCellAsString(i, 0), t.getCellAsString(i, 7));
        }

        try {
            BufferedReader br = IOUtils.getTextReader(delfileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Taxa\tSub\tRatio\tGroup\tSubspecies\tGroupID");
            bw.newLine();
            String temp = br.readLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) { //读的是有害突变那一列
                l = PStringUtils.fastSplit(temp);
                cnt++;
                String taxa = l.get(0);
                String sub = l.get(1);
                String ratiodel = l.get(7);
                Record r = new Record(taxa,sub);
               int index = Collections.binarySearch(taxasubl,r);
               if(index > -1){
                   double ratiosyn = Double.parseDouble(hmratio.get(cnt));
                   double rr = (double)Double.parseDouble(ratiodel)/ratiosyn;
                   bw.write(taxa + "\t" + sub + "\t" + String.format("%.3f",rr)+"\t"+taxaGroupMap.get(taxa)+"\t"+taxaSubMap.get(taxa)+"\t"+taxaGroupIDMap.get(taxa));
                   bw.newLine();
               }
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

    /**
     * 在taxaList文件中添加index，便于查看最详细的亚群（如 云南小麦 西藏半野生小麦的 mutation burden）
     */
    public void mkIndexforBurden(){
        try {
            String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/000_prepareData/001_input/taxaList.txt";
            String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/000_prepareData/001_input/taxaList_addIndexBurden.txt";
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = br.readLine();
            bw.write(temp + "\tIndexforMutationBurden");
            bw.newLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            new AoFile().readheader(infileS);
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                String in = l.get(6);
                String goal = in.split("_")[0];
                bw.write(temp + "\t" + goal);
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
    //根据最终生成的文件，进行 A B D sub的合并
    public void mergeFinalfilebySub(){
        //change
//        String infileS = "";

//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/003_delCount_highDepth/001_/additiveDeleterious_vmap2_highDepth_bychr.txt";
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/003_delCount_highDepth/002_/additiveDeleterious_vmap2_highDepth_bychr.txt";
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/003_delCount_highDepth/003_/additiveDeleterious_vmap2_highDepth_bychr.txt";
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/003_delCount_highDepth/000_/additiveDeleterious_vmap2_highDepth_bychr.txt";

        String splitDirS = new File(infileS).getParent()+"/split";
        new File(splitDirS).mkdirs();

        String outfileS = infileS.replaceFirst("_bychr.txt", "_bysub.txt");

        //查看有多少个taxa
        RowTable<String> t = new RowTable<>(infileS);
        Set<String> s = new HashSet<>(t.getColumn(0));
        String[] taxa = s.toArray(new String[s.size()]);
        Arrays.sort(taxa);
        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            String header = br.readLine();
            BufferedWriter[] bw = new BufferedWriter[taxa.length];
            for (int i = 0; i < taxa.length; i++) {
                String tempoutS = new File(splitDirS,taxa[i]+".txt").getAbsolutePath();
                bw[i] = IOUtils.getTextWriter(tempoutS);
                bw[i].write(header);
                bw[i].newLine();
            }

            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String ta = l.get(0);
                int index = Arrays.binarySearch(taxa,ta);
                bw[index].write(temp);
                bw[index].newLine();
                cnt++;

            }
            br.close();
            for (int i = 0; i < taxa.length; i++) {
                bw[i].flush();
                bw[i].close();
            }
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }


        //不变的
        String taxaGroupFileS = "/Users/Aoyue/project/wheatVMapII/001_germplasm/GermplasmDB/001_toFeiLu/004_wheatVMapII_germplasmInfoAddSubspecies_20191225.txt";
        String taxaGroup2FileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/000_prepareData/001_input/taxaList.txt";
        t = new RowTable (taxaGroupFileS);
        HashMap<String, String> taxaGroupMap = new HashMap();
        HashMap<String, String> taxaSubMap = new HashMap();
        HashMap<String, String> taxaGroupIDMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            taxaGroupMap.put(t.getCellAsString(i, 4), t.getCellAsString(i, 26));
            taxaSubMap.put(t.getCellAsString(i, 4), t.getCellAsString(i, 27));
        }
        t = new RowTable (taxaGroup2FileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            taxaGroupIDMap.put(t.getCellAsString(i, 0), t.getCellAsString(i, 7));
        }

        HashMap<String,Integer> h = new HashMap<>();
        h.put("A",0);
        h.put("B",1);
        h.put("D",2);

        String outDirS = null;
        for (int i = 0; i < taxa.length; i++) {
            try {
                outDirS = new File(splitDirS).getParent() + "/split_merge";
                new File(outDirS).mkdirs();

                String inS = new File(splitDirS,taxa[i]+".txt").getAbsolutePath();
                String outS = new File(outDirS,taxa[i]+"_bysub.txt").getAbsolutePath();
                BufferedReader br = IOUtils.getTextReader(inS);
                BufferedWriter bw = IOUtils.getTextWriter(outS);
                bw.write(br.readLine().replaceFirst("Chr","Sub"));
                bw.newLine();
                TDoubleArrayList[] derivedDelList = new TDoubleArrayList[3];
                TDoubleArrayList[] genotypeList = new TDoubleArrayList[3];
                for (int j = 0; j < derivedDelList.length; j++) {
                    derivedDelList[j] = new TDoubleArrayList();
                    genotypeList[j] = new TDoubleArrayList();

                }
                String temp = null;
                List<String> l = new ArrayList<>();
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    cnt++;
                    int chr = Integer.parseInt(l.get(1));
                    double d = Double.parseDouble(l.get(2));
                    double k = Double.parseDouble(l.get(3));
                    String sub = this.convertChrtoSub(chr);
                    derivedDelList[h.get(sub)].add(d);
                    genotypeList[h.get(sub)].add(k);
                }
                br.close();


                //求和
                DescriptiveStatistics[] d = new DescriptiveStatistics[3];
                DescriptiveStatistics[] dd = new DescriptiveStatistics[3];

                Double[] ratio = new Double[3];
                for (int j = 0; j < derivedDelList.length; j++) {
                    d[j] = new DescriptiveStatistics(derivedDelList[j].toArray());
                    dd[j] = new DescriptiveStatistics(genotypeList[j].toArray());
                    ratio[j] = d[j].getSum()/dd[j].getSum();

                }

                HashMap<Integer,String> hhh = new HashMap<>();
                hhh.put(0,"A");
                hhh.put(1,"B");
                hhh.put(2,"D");


                for (int j = 0; j < derivedDelList.length; j++) {
                    if(dd[j].getSum()==0)continue;
                    bw.write(taxa[i] + "\t" + hhh.get(j) + "\t" + String.format("%.1f",d[j].getSum()) + "\t" + String.format("%.0f",dd[j].getSum())
                            + "\t" + taxaGroupMap.get(taxa[i]) + "\t" + taxaSubMap.get(taxa[i]) + "\t" + taxaGroupIDMap.get(taxa[i]) + "\t" + String.format("%.3f",ratio[j]));
                    bw.newLine();
                }

                bw.flush();
                bw.close();

            } catch (Exception e) {
                e.printStackTrace();
            }

        }

        new CountSites().mergeTxtbysuffix(outDirS,outfileS,".txt");
        new File(splitDirS).delete();
        new File(outDirS).delete();

    }


    /**
     * 根据染色体号，返回所在亚基因组
     * @param a
     * @return
     */
    public String convertChrtoSub(int a){
        int[] arra = {1, 2, 7, 8, 13, 14, 19, 20, 25, 26, 31, 32, 37, 38};
        int[] arrb = {3, 4, 9, 10, 15, 16, 21, 22, 27, 28, 33, 34, 39, 40};
        int[] arrd = {5, 6, 11, 12, 17, 18, 23, 24, 29, 30, 35, 36, 41, 42};
        HashMap<Integer, String> hml = new HashMap<>();
        Arrays.sort(arra);
        Arrays.sort(arrb);
        Arrays.sort(arrd);
        for (int i = 0; i < arra.length; i++) {
            hml.put(arra[i], "A");
            hml.put(arrb[i], "B");
            hml.put(arrd[i], "D");
        }
        String subgenome = hml.get(a);
        return subgenome;
    }


    private void countDeleteriousVMapIIHighDepth () {
        //不变的文件
        String taxaSummaryFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/002_VMapIIDepth/taxaDepth_summary.txt";
        String taxaGroupFileS = "/Users/Aoyue/project/wheatVMapII/001_germplasm/GermplasmDB/001_toFeiLu/004_wheatVMapII_germplasmInfoAddSubspecies_20191225.txt";
        String taxaGroup2FileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/000_prepareData/001_input/taxaList.txt";

        //可变的文件
//        String addInfileS = "";
//        String addOutfileS = "";

//        String addInfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/002_/additiveDeleterious_vmap2_bychr.txt";
//        String addOutfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/003_delCount_highDepth/002_/additiveDeleterious_vmap2_highDepth_bychr.txt";

//        String addInfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/003_/additiveDeleterious_vmap2_bychr.txt";
//        String addOutfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/003_delCount_highDepth/003_/additiveDeleterious_vmap2_highDepth_bychr.txt";

        String addInfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/000_/additiveDeleterious_vmap2_bychr.txt";
        String addOutfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/003_delCount_highDepth/000_/additiveDeleterious_vmap2_highDepth_bychr.txt";



//        String addInfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/additiveDeleterious_vmap2.txt";
//        String addInfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/001_/additiveDeleterious_vmap2_bychr.txt";
//        String addOutfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/003_delCount_highDepth/001_/additiveDeleterious_vmap2_highDepth_bychr.txt";
//        String recInfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/001_/recessiveDeleterious_vmap2_bychr.txt";
//        String recOutfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/003_delCount_highDepth/001_/recessiveDeleterious_vmap2_highDepth_bychr.txt";
//        String taxaGroupFileS = "/Users/Aoyue/project/wheatVMapII/001_germplasm/GermplasmDB/001_toFeiLu/004_wheatVMapII_germplasmInfoAddSubspecies_20191225.txt";
//        String taxaGroup2FileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/000_prepareData/001_input/taxaList.txt";


//        String taxaSummaryFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/002_VMapIIDepth/taxaDepth_summary.txt";
//        String addInfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/additiveSyn_vmap2.txt";
//        String addOutfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/003_delCount_highDepth/additiveSyn_vmap2_highDepth.txt";
//        String recInfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/recessiveSyn_vmap2.txt";
//        String recOutfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/003_delCount_highDepth/recessiveSyn_vmap2_highDepth.txt";
//        String taxaGroupFileS = "/Users/Aoyue/project/wheatVMapII/001_germplasm/GermplasmDB/001_toFeiLu/004_wheatVMapII_germplasmInfoAddSubspecies_20191225.txt";
//        String taxaGroup2FileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/000_prepareData/001_input/taxaList.txt";



        new AoFile().readheader(taxaGroupFileS);
        double depthCut = 3;
        RowTable t = new RowTable (taxaSummaryFileS);
        ArrayList<String> taxaList = new ArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (t.getCellAsDouble(i, 2) < depthCut) continue;
            taxaList.add(t.getCellAsString(i, 0));
        }
        String[] taxa = taxaList.toArray(new String[taxaList.size()]); //taxa是具有高深度的taxa列表
        Arrays.sort(taxa);
        t = new RowTable (taxaGroupFileS);
        HashMap<String, String> taxaGroupMap = new HashMap();
        HashMap<String, String> taxaSubMap = new HashMap();
        HashMap<String, String> taxaGroupIDMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            taxaGroupMap.put(t.getCellAsString(i, 4), t.getCellAsString(i, 26));
            taxaSubMap.put(t.getCellAsString(i, 4), t.getCellAsString(i, 27));
        }

        t = new RowTable (taxaGroup2FileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            taxaGroupIDMap.put(t.getCellAsString(i, 0), t.getCellAsString(i, 7));
        }

        t = new RowTable (addInfileS);
        List<String> headerl = t.getHeader();

        try {
            BufferedWriter bw = IOUtils.getTextWriter(addOutfileS);
            for (int i = 0; i < headerl.size(); i++) {
                bw.write(headerl.get(i) + "\t");
            }
            bw.write("Group\tSubspecies\tGroupID\tRatio");
//            bw.write("Taxa\tDeleteriousCountPerHaplotype\tSiteCountWithMinDepth\tGroup\tSubspecies\tGroupID\tRatio");
            bw.newLine();
            for (int i = 0; i < t.getRowNumber(); i++) { //这里的t是 addInfileS
                int index = Arrays.binarySearch(taxa, t.getCellAsString(i, 0));
                if (index < 0) continue;
                if(taxaGroupMap.get(taxa[index]).equals("ExclusionHexaploid") || taxaGroupMap.get(taxa[index]).equals("ExclusionTetraploid")) continue;
                double genotypesite = Double.valueOf(t.getCellAsDouble(i, 2));
                if(genotypesite == 0) continue;
                double ratio = Double.valueOf(t.getCellAsDouble(i, 2))/Double.valueOf(t.getCellAsDouble(i, 3));
                for (int j = 0; j < t.getColumnNumber(); j++) {
                    bw.write(t.getCellAsString(i,j) + "\t");
                }
                bw.write(taxaGroupMap.get(taxa[index])
                        +"\t"+taxaSubMap.get(taxa[index])
                        +"\t"+taxaGroupIDMap.get(taxa[index])
                        +"\t"+String.format("%.3f",ratio));
                bw.newLine();
            }
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }


    private void countDeleteriousVMapIIHighDepth_old () {
        String taxaSummaryFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/002_VMapIIDepth/taxaDepth_summary.txt";
        String addInfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/001_/additiveDeleterious_vmap2.txt";
        String addOutfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/003_delCount_highDepth/001_/additiveDeleterious_vmap2_highDepth.txt";
        String recInfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/001_/recessiveDeleterious_vmap2.txt";
        String recOutfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/003_delCount_highDepth/001_/recessiveDeleterious_vmap2_highDepth.txt";
        String taxaGroupFileS = "/Users/Aoyue/project/wheatVMapII/001_germplasm/GermplasmDB/001_toFeiLu/004_wheatVMapII_germplasmInfoAddSubspecies_20191225.txt";
        String taxaGroup2FileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/000_prepareData/001_input/taxaList.txt";
        //1.54=3x, 2.75 = 5x

//        String taxaSummaryFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/002_VMapIIDepth/taxaDepth_summary.txt";
//        String addInfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/additiveSyn_vmap2.txt";
//        String addOutfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/003_delCount_highDepth/additiveSyn_vmap2_highDepth.txt";
//        String recInfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/recessiveSyn_vmap2.txt";
//        String recOutfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/003_delCount_highDepth/recessiveSyn_vmap2_highDepth.txt";
//        String taxaGroupFileS = "/Users/Aoyue/project/wheatVMapII/001_germplasm/GermplasmDB/001_toFeiLu/004_wheatVMapII_germplasmInfoAddSubspecies_20191225.txt";
//        String taxaGroup2FileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/000_prepareData/001_input/taxaList.txt";
        //1.54=3x, 2.75 = 5x
        new AoFile().readheader(taxaGroupFileS);
        double depthCut = 3;
        RowTable t = new RowTable (taxaSummaryFileS);
        ArrayList<String> taxaList = new ArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (t.getCellAsDouble(i, 2) < depthCut) continue;
            taxaList.add(t.getCellAsString(i, 0));
        }
        String[] taxa = taxaList.toArray(new String[taxaList.size()]);
        Arrays.sort(taxa);
        t = new RowTable (taxaGroupFileS);
        HashMap<String, String> taxaGroupMap = new HashMap();
        HashMap<String, String> taxaSubMap = new HashMap();
        HashMap<String, String> taxaGroupIDMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            taxaGroupMap.put(t.getCellAsString(i, 4), t.getCellAsString(i, 26));
            taxaSubMap.put(t.getCellAsString(i, 4), t.getCellAsString(i, 27));
        }

        t = new RowTable (taxaGroup2FileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            taxaGroupIDMap.put(t.getCellAsString(i, 0), t.getCellAsString(i, 7));
        }

        t = new RowTable (addInfileS);
        try {
            BufferedWriter bw = IOUtils.getTextWriter(addOutfileS);
            bw.write("Taxa\tDeleteriousCountPerHaplotype\tSiteCountWithMinDepth\tGroup\tSubspecies\tGroupID\tRatio");
            bw.newLine();
            for (int i = 0; i < t.getRowNumber(); i++) {
                int index = Arrays.binarySearch(taxa, t.getCellAsString(i, 0));
                if (index < 0) continue;
                if(taxaGroupMap.get(taxa[index])=="ExclusionHexaploid" || taxaGroupMap.get(taxa[index])=="ExclusionTetraploid") continue;
                double ratio = Double.valueOf(t.getCellAsDouble(i, 1))/Double.valueOf(t.getCellAsDouble(i, 2));
                bw.write(taxa[index]+"\t"+t.getCellAsString(i, 1)+"\t"+t.getCellAsString(i, 2)+"\t"+taxaGroupMap.get(taxa[index])
                        +"\t"+taxaSubMap.get(taxa[index])+"\t"+taxaGroupIDMap.get(taxa[index])
                        +"\t"+String.format("%.3f",ratio));
                bw.newLine();
            }
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        t = new RowTable (recInfileS);
        try {
            BufferedWriter bw = IOUtils.getTextWriter(recOutfileS);
            bw.write("Taxa\tDeleteriousCountPerLine\tSiteCountWithMinDepth\tGroup\tRatio");
            bw.newLine();
            for (int i = 0; i < t.getRowNumber(); i++) {
                int index = Arrays.binarySearch(taxa, t.getCellAsString(i, 0));
                if (index < 0) continue;
                double ratio = Double.valueOf(t.getCellAsDouble(i, 1))/Double.valueOf(t.getCellAsDouble(i, 2));
                bw.write(taxa[index]+"\t"+t.getCellAsString(i, 1)+"\t"+t.getCellAsString(i, 2)+"\t"+taxaGroupMap.get(taxa[index])
                        +"\t"+taxaSubMap.get(taxa[index])
                        +"\t"+String.format("%.3f",ratio));
                bw.newLine();
            }
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void mergeDepthSummary(){
        String infileS =  "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/002_VMapIIDepth/taxaDepth_Asub.summary.txt";
        String infile2S = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/002_VMapIIDepth/taxaDepth_Dsub.summary.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/002_VMapIIDepth/taxaDepth_summary.txt";
        HashMap<String,String> hm1 = new AoFile().getHashMapStringKey(infileS,0,2);
        HashMap<String,String> hm2 = new AoFile().getHashMapStringKey(infile2S,0,2);
        Set<String> s = new HashSet<>();
        for (int i = 0; i < hm1.size(); i++) {
            s.addAll(hm1.keySet());
        }
        for (int i = 0; i < hm2.size(); i++) {
            s.addAll(hm2.keySet());
        }
        System.out.println(s.size());
        List<String> l = new ArrayList<>(s);
        Collections.sort(l);

        List<String> l1 = new ArrayList<>(hm1.keySet());
        List<String> l2 = new ArrayList<>(hm2.keySet());
        Collections.sort(l1);
        Collections.sort(l2);

        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Taxa\tID\tMeanDepth");
            bw.newLine();
            String temp = null;
            int cnt = 0;
            String value = null;
            for (int i = 0; i < l.size(); i++) {
                int index1 = Collections.binarySearch(l1,l.get(i));
                int index2 = Collections.binarySearch(l2,l.get(i));
                if(index1 > -1 && (index2 > -1)){ //2个文件都存在，就写第一个文件
                    bw.write(l.get(i) + "\t" + (i+1) + "\t" + hm1.get(l.get(i)));
                    bw.newLine();
                }
                if(index1 > -1 && (index2 < 0)){
                    bw.write(l.get(i) + "\t" + (i+1) + "\t" + hm1.get(l.get(i)));
                    bw.newLine();
                }
                if(index1 < 0 && (index2 > -1)){
                    bw.write(l.get(i) + "\t" + (i+1) + "\t" + hm2.get(l.get(i)));
                    bw.newLine();
                }
            }

            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public void mkDepthSummary () {
        //思想：对每一个taxa做统计，每读进一个taxa，就统计所有深度的和，再除以表格的行数，也即是抽查的位点数，最终得出总得深度除以位点数，得出平均每个位点的深度。
//        String taxaDepthDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/002_VMapIIDepth/taxa_Asub";
//        String taxaSummaryFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/002_VMapIIDepth/taxaDepth_Asub.summary.txt";

        String taxaDepthDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/002_VMapIIDepth/taxa_Dsub";
        String taxaSummaryFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/002_VMapIIDepth/taxaDepth_Dsub.summary.txt";


        File[] fs = new File (taxaDepthDirS).listFiles();
        Arrays.sort(fs);
        try {
            BufferedWriter bw = IOUtils.getTextWriter(taxaSummaryFileS);
            bw.write("Taxa\tID\tMeanDepth");
            bw.newLine();
            for (int i = 0; i < fs.length; i++) {
                RowTable t = new RowTable (fs[i].getAbsolutePath());
                String taxaName = String.valueOf(t.getHeader().get(0)).replaceFirst("_siteDepth", "");
                double value = 0;
                for (int j = 0; j < t.getRowNumber(); j++) {
                    value+=t.getCellAsDouble(j, 0);
                }
                double dd = (double)value/t.getRowNumber();
                bw.write(taxaName+"\t"+String.valueOf(i+1)+"\t"+String.format("%.2f", dd));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void mkDepthOfVMapII(){
//        String vcfFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/024_subsetVCF_maf0.01byPop/002_merged/chr.Asub.maf0.01byPop.vcf.gz";
//        String hmpInfoFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/014_merge/chr_A.SNP_anno.txt.gz";
//        String taxaDepthDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/002_VMapIIDepth/taxa_Asub";

        String vcfFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/024_subsetVCF_maf0.01byPop/002_merged/chr.Dsub.maf0.01byPop.vcf.gz";
        String hmpInfoFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/014_merge/chr_D.SNP_anno.txt.gz";
        String taxaDepthDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/002_VMapIIDepth/taxa_Dsub";

        int snpNum = 0;
        int size = 500000;
        try {
            BufferedReader br = IOUtils.getTextGzipReader(hmpInfoFileS);
            String temp = br.readLine();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                cnt++;
            }
            snpNum = cnt;
            int[] indices = new int[size];
            for (int i = 0; i < size; i++) {
                indices[i] = (int)(Math.random()*snpNum);
            }
            Arrays.sort(indices);
            br = IOUtils.getTextGzipReader(vcfFileS);
            while ((temp = br.readLine()).startsWith("##")) {}
            List<String> l = PStringUtils.fastSplit(temp, "\t");
            String[] taxa = new String[l.size()-9];
            for (int i = 0; i < taxa.length; i++) {
                taxa[i] = l.get(i+9);
            }
            TIntArrayList[] depthList = new TIntArrayList[taxa.length];
            for (int i = 0; i < taxa.length; i++) depthList[i] = new TIntArrayList();
            cnt = 0;
            while ((temp = br.readLine()) != null) {
                cnt++; // 对snp开始计数
                if (cnt%1000000 == 0) System.out.println(String.valueOf(cnt)+" lines");
                int idx = Arrays.binarySearch(indices, cnt-1);
                if (idx < 0) continue;
                l = PStringUtils.fastSplit(temp, "\t");
                for (int i = 0; i < taxa.length; i++) {
                    String genoS = l.get(i+9);
                    if (genoS.startsWith(".")) {
                        depthList[i].add(0);
                        continue;
                    }
                    List<String> ll = PStringUtils.fastSplit(genoS, ":");
                    List<String> lll = PStringUtils.fastSplit(ll.get(1), ",");
                    int depth = Integer.valueOf(lll.get(0))+Integer.valueOf(lll.get(1));
                    depthList[i].add(depth);
                }
            }
            for (int i = 0; i < taxa.length; i++) {
                String outfileS = new File (taxaDepthDirS, PStringUtils.getNDigitNumber(3, i+1)+"depth.txt").getAbsolutePath();
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                bw.write(taxa[i]+"_siteDepth");
                bw.newLine();
                int[] depth = depthList[i].toArray();
                for (int j = 0; j < depth.length; j++) {
                    bw.write(String.valueOf(depth[j]));
                    bw.newLine();
                }
                bw.flush();
                bw.close();
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }

    }

    /**
     * 新的数据库分析，添加了dirivedSIFT的列信息,并且将每条染色体的情况都统计下来
     */
    public void countDeleteriousVMapII_byChr() {

//        String delVCFDirS = ""; //有害变异的VCF文件路径
//        String deleFileS = ""; //有害变异信息库
//        String addCountFileS = ""; //有害变异加性模型输出文件
//        String recCountFileS = ""; //有害变异隐形模型输出文件

//        String delVCFDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/002_vcf/001_/del"; //有害变异的VCF文件路径
//        String deleFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/001_snp/001_/del_merged/chr_SNP_anno_del.txt.gz"; //有害变异信息库
//        String addCountFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/001_/additiveDeleterious_vmap2_bychr.txt"; //有害变异加性模型输出文件
//        String recCountFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/001_/recessiveDeleterious_vmap2_bychr.txt"; //有害变异隐形模型输出文件

//        String delVCFDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/002_vcf/002_/del"; //有害变异的VCF文件路径
//        String deleFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/001_snp/002_/del_merged/chr_SNP_anno_del.txt.gz"; //有害变异信息库
//        String addCountFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/002_/additiveDeleterious_vmap2_bychr.txt"; //有害变异加性模型输出文件
//        String recCountFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/002_/recessiveDeleterious_vmap2_bychr.txt"; //有害变异隐形模型输出文件

//        String delVCFDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/002_vcf/003_/syn"; //有害变异的VCF文件路径
//        String deleFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/001_snp/003_/syn_merged/chr_SNP_anno_syn.txt.gz"; //有害变异信息库
//        String addCountFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/003_/additiveDeleterious_vmap2_bychr.txt"; //有害变异加性模型输出文件
//        String recCountFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/003_/recessiveDeleterious_vmap2_bychr.txt"; //有害变异隐形模型输出文件

        String delVCFDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/002_vcf/000_/del"; //有害变异的VCF文件路径
        String deleFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/001_snp/000_/del_merged/chr_SNP_anno_del.txt.gz"; //有害变异信息库
        String addCountFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/000_/additiveDeleterious_vmap2_bychr.txt"; //有害变异加性模型输出文件
        String recCountFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/000_/recessiveDeleterious_vmap2_bychr.txt"; //有害变异隐形模型输出文件


        int minDepth = 2;//inclusive
        int chrNum = 42;
        ArrayList<Integer> chrList = new ArrayList();
        for (int i = 0; i < chrNum; i++) {
            chrList.add(i + 1);
        }


        RowTable<String> t = new RowTable(deleFileS);
        TIntArrayList[] posList = new TIntArrayList[chrNum];
        TCharArrayList[] charList = new TCharArrayList[chrNum];
        for (int i = 0; i < chrNum; i++) { //集合类数组，要初始化每一个list
            posList[i] = new TIntArrayList();
            charList[i] = new TCharArrayList();
        }
        new AoFile().readheader(deleFileS);
        String derivedAllele = null;
        for (int i = 0; i < t.getRowNumber(); i++) {
            int index = t.getCellAsInteger(i, 1) - 1; //染色体号的索引
            String ancestralAllele = t.getCell(i, 15); //不同的数据库，这一列的信息不一样，千万要注意
            String majorAllele = t.getCell(i, 5);
            String minorAllele = t.getCell(i, 6);
            if (ancestralAllele.equals(majorAllele)) {
                derivedAllele = minorAllele;
                posList[index].add(t.getCellAsInteger(i, 2)); //将包含有derived allele的位点添加到Poslist
                charList[index].add(derivedAllele.charAt(0)); //返回的是char类型的字符
            }
            if (ancestralAllele.equals(minorAllele)) {
                derivedAllele = majorAllele;
                posList[index].add(t.getCellAsInteger(i, 2));
                charList[index].add(derivedAllele.charAt(0));
            }
        }
        int[][] delePos = new int[chrNum][];
        char[][] deleChar = new char[chrNum][];
        for (int i = 0; i < chrNum; i++) { //将每一个list转化为数组
            delePos[i] = posList[i].toArray();
            deleChar[i] = charList[i].toArray();
            Arrays.sort(delePos[i]);
        }


        String vmap2TaxaList = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/004_vmap2.1_TaxaList/TaxaList.txt";
        List<String> taxaList = new AoFile().getStringListwithoutHeader(vmap2TaxaList, 0);
        String[] taxa = taxaList.toArray(new String[taxaList.size()]);
        Arrays.sort(taxa);


        int taxaNum = taxa.length;
        double[][] addCount = new double[chrNum][taxa.length];
        int[][] recCount = new int[chrNum][taxa.length];
        int[][] siteWithMinDepthCount = new int[chrNum][taxa.length]; //每个taxa的深度

        chrList.parallelStream().forEach(chr -> {
            String delVmapFileS = "chr" + PStringUtils.getNDigitNumber(3, chr) + "_del_vmap2.1.vcf.gz";
//            String delVmapFileS = "chr" + PStringUtils.getNDigitNumber(3, chr) + "_syn_vmap2.1.vcf.gz"; /////////此处需要修改

            delVmapFileS = new File(delVCFDirS, delVmapFileS).getAbsolutePath();
            BufferedReader br = IOUtils.getTextGzipReader(delVmapFileS);
            int chrIndex = chr - 1;
            try {
                String temp = null;
                int cnt = 0;
                HashMap<String, Integer> hmtaxainVCFindex = new HashMap<>();
                HashMap<String, Integer> hmtaxainTaxaindex = new HashMap<>();
                //这里涉及2个数组概念， taxa and taxainVCFfile,即总的taxa和在VCF文件中的taxa
                String[] taxainVCFfile = null;
                List<String> taxainVCFlist = new ArrayList();
                while ((temp = br.readLine()) != null) {

                    if (temp.startsWith("##")) {
                        continue;
                    }
                    if (temp.startsWith("#CHROM")) {//说明进入taxa信息列
                        List<String> l = PStringUtils.fastSplit(temp, "\t");
                        taxainVCFfile = new String[l.size() - 9];
                        for (int i = 9; i < l.size(); i++) {
                            taxainVCFfile[i - 9] = l.get(i);
                            hmtaxainVCFindex.put(taxainVCFfile[i - 9], i); //第0个taxa的genotype在第9行，第一个taxa的genotype在第10行，依次类推；
                        }
                        // ************** 在总的taxa中，搜索VCF文件中的taxa的index
                        for (int i = 0; i < taxainVCFfile.length; i++) { //找到taxa
                            int index = Arrays.binarySearch(taxa, taxainVCFfile[i]);
                            if (index > -1) {
                                hmtaxainTaxaindex.put(taxainVCFfile[i], index);
                            }
                        }
                        taxainVCFlist = Arrays.asList(taxainVCFfile);
                    }
                    if (!temp.startsWith("#")) {
                        cnt++;
                        if (cnt % 1000 == 0) {
                            System.out.println(String.valueOf(cnt) + " lines on chr " + String.valueOf(chr));
                        }

                        List<String> l = PStringUtils.fastSplit(temp.substring(0, 100), "\t");
                        int pos = Integer.valueOf(l.get(1));
                        int index = Arrays.binarySearch(delePos[chrIndex], pos);
                        if (index < 0) {
                            continue; //前面建立的 delePos[][] 和deleChar[][] 都是为现在在vcf文件中找位置贡献的，不是有害突变的位点，都过滤。*************************************************************
                        }
                        l = PStringUtils.fastSplit(temp, "\t");
                        int[] idx = new int[2];
                        if (l.get(3).charAt(0) == deleChar[chrIndex][index]) { //如果ref allele = deleChar allele
                            idx[0] = 0; //等于0则代表是 derived allele
                            idx[1] = 1;
                        } else {
                            idx[0] = 1;
                            idx[1] = 0;
                        }


                        //合计642个taxa，在A Bgenome中只有606（419+187）个，在Dgenome中只有455（419+36）个，我们要找到每个VCF文件中的genotype所对应的taxa的index
                        //code:
                        for (int i = 0; i < taxainVCFlist.size(); i++) {
                            int genotypeIndex = hmtaxainVCFindex.get(taxainVCFlist.get(i)); //获取该taxa的基因型所在的列的索引
                            int taxaIndex = hmtaxainTaxaindex.get(taxainVCFlist.get(i)); //获取该taxa所在总的642个数组中的索引，为后续写文件进行统计
                            String genoS = l.get(genotypeIndex); //指的是 GT:AD:GL 信息
                            if (genoS.startsWith(".")) {
                                continue; //如果以.开头，说明没有基因型信息，此位点没有测到。
                            }
                            List<String> ll = PStringUtils.fastSplit(genoS, ":"); //分开为GT   AD   GL三类

                            List<String> lll = PStringUtils.fastSplit(ll.get(1), ","); //lll指将AD提取出来，并以"，"号分割。如 0/0:1,0:0,3,28中 ，1，0分别代表ref和alt的测序深度
                            int depth = Integer.valueOf(lll.get(0)) + Integer.valueOf(lll.get(1)); //总得测序深度等于 ref + alt
                            if (depth < minDepth) {
                                continue; //最小测序深度是2，如果小于2，则弃用
                            }
                            lll = PStringUtils.fastSplit(ll.get(0), "/"); //这里lll指的是基因型GT，lll被重新赋值，之前代表的是AD
                            int v1 = Integer.valueOf(lll.get(0)); //v1等于 ref
                            int v2 = Integer.valueOf(lll.get(1)); // v2 等于 alt
                            int sum = 0;
                            //如果ref是derived allele，那么idx[0]=0;idx[1]=1. 当时0/0时，sum=2; 0/1时，sum=1; 1/1时， sum=0.
                            //如果alt是derived allele，那么idx[0]=1;idx[1]=0. 当时0/0时，sum=0; 0/1时，sum=1; 1/1时， sum=2.
                            //
                            //
                            if (v1 == idx[0]) { //
                                sum++;
                            }
                            if (v2 == idx[0]) {
                                sum++;
                            }
                            if (sum == 0) {
                            } else if (sum == 1) {
                                addCount[chrIndex][taxaIndex] += 0.5;
                            } else {
                                addCount[chrIndex][taxaIndex] += 1;
                                recCount[chrIndex][taxaIndex] += 1;
                            }
                            siteWithMinDepthCount[chrIndex][taxaIndex]++; //taxa有多少个有害突变位点
                        }
                    }
                }
                br.close();
            } catch (Exception e) {
                e.printStackTrace();
            }

            /**
             * 进行每个样品每条染色体的计算
             */
        });

        try {
            BufferedWriter bw = IOUtils.getTextWriter(addCountFileS);
            bw.write("Taxa\tChr\tDeleteriousCountPerHaplotype\tSiteCountWithMinDepth"); //每个taxa有多少个加性效应的derivedAllele 每个taxa在del库中含有基因型的个数
            bw.newLine();
            for (int i = 0; i < addCount.length; i++) {
                String chr = String.valueOf(i + 1);
                for (int j = 0; j < addCount[0].length; j++) {
                    bw.write(taxa[j] + "\t" + chr + "\t" + String.valueOf(addCount[i][j]) + "\t" + String.valueOf(siteWithMinDepthCount[i][j]));
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();

            bw = IOUtils.getTextWriter(recCountFileS);
            bw.write("Taxa\tChr\tDeleteriousCountPerHaplotype\tSiteCountWithMinDepth"); //每个taxa有多少个加性效应的derivedAllele 每个taxa在del库中含有基因型的个数
            bw.newLine();
            for (int i = 0; i < addCount.length; i++) {
                String chr = String.valueOf(i + 1);
                for (int j = 0; j < addCount[0].length; j++) {
                    bw.write(taxa[j] + "\t" + chr + "\t" + String.valueOf(recCount[i][j]) + "\t" + String.valueOf(siteWithMinDepthCount[i][j]));
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();

        } catch (Exception e) {
            e.printStackTrace();
        }
    }


    /**
     * 新的数据库分析，添加了dirivedSIFT的列信息
     */
    public void countDeleteriousVMapII() {
//        String delVCFDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/002_vcf/del"; //有害变异的VCF文件路径
//        String deleFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/001_snp/del_merged/chr_SNP_anno_del.txt.gz"; //有害变异信息库
//        String addCountFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/additiveDeleterious_vmap2.txt"; //有害变异加性模型输出文件
//        String recCountFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/recessiveDeleterious_vmap2.txt"; //有害变异隐形模型输出文件
//
        String delVCFDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/002_vcf/syn"; //有害变异的VCF文件路径
        String deleFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/001_snp/syn_merged/chr_SNP_anno_syn.txt.gz"; //有害变异信息库
        String addCountFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/additiveSyn_vmap2.txt"; //有害变异加性模型输出文件
        String recCountFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/recessiveSyn_vmap2.txt"; //有害变异隐形模型输出文件


        int minDepth = 2;//inclusive
        int chrNum = 42;
        ArrayList<Integer> chrList = new ArrayList();
        for (int i = 0; i < chrNum; i++) {
            chrList.add(i + 1);
        }


        RowTable<String> t = new RowTable(deleFileS);
        TIntArrayList[] posList = new TIntArrayList[chrNum];
        TCharArrayList[] charList = new TCharArrayList[chrNum];
        for (int i = 0; i < chrNum; i++) { //集合类数组，要初始化每一个list
            posList[i] = new TIntArrayList();
            charList[i] = new TCharArrayList();
        }
        new AoFile().readheader(deleFileS);
        String derivedAllele = null;
        for (int i = 0; i < t.getRowNumber(); i++) {
            int index = t.getCellAsInteger(i, 1) - 1; //染色体号的索引
            String ancestralAllele = t.getCell(i, 15);
            String majorAllele = t.getCell(i, 5);
            String minorAllele = t.getCell(i, 6);
            if (ancestralAllele.equals(majorAllele)) {
                derivedAllele = minorAllele;
                posList[index].add(t.getCellAsInteger(i, 2)); //将包含有derived allele的位点添加到Poslist
                charList[index].add(derivedAllele.charAt(0)); //返回的是char类型的字符
            }
            if (ancestralAllele.equals(minorAllele)) {
                derivedAllele = majorAllele;
                posList[index].add(t.getCellAsInteger(i, 2));
                charList[index].add(derivedAllele.charAt(0));
            }
        }
        int[][] delePos = new int[chrNum][];
        char[][] deleChar = new char[chrNum][];
        for (int i = 0; i < chrNum; i++) { //将每一个list转化为数组
            delePos[i] = posList[i].toArray();
            deleChar[i] = charList[i].toArray();
            Arrays.sort(delePos[i]);
        }


        String vmap2TaxaList = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/004_vmap2.1_TaxaList/TaxaList.txt";
        List<String> taxaList = new AoFile().getStringListwithoutHeader(vmap2TaxaList, 0);
        String[] taxa = taxaList.toArray(new String[taxaList.size()]);
        Arrays.sort(taxa);


        int taxaNum = taxa.length;
        double[] addCount = new double[taxa.length];
        int[] recCount = new int[taxa.length];
        int[] siteWithMinDepthCount = new int[taxa.length]; //每个taxa的深度

        chrList.parallelStream().forEach(chr -> {
//            String delVmapFileS = "chr" + PStringUtils.getNDigitNumber(3, chr) + "_del_vmap2.1.vcf.gz";
            String delVmapFileS = "chr" + PStringUtils.getNDigitNumber(3, chr) + "_syn_vmap2.1.vcf.gz"; /////////此处需要修改

            delVmapFileS = new File(delVCFDirS, delVmapFileS).getAbsolutePath();
            BufferedReader br = IOUtils.getTextGzipReader(delVmapFileS);
            int chrIndex = chr - 1;
            try {
                String temp = null;
                int cnt = 0;
                HashMap<String, Integer> hmtaxainVCFindex = new HashMap<>();
                HashMap<String, Integer> hmtaxainTaxaindex = new HashMap<>();
                //这里涉及2个数组概念， taxa and taxainVCFfile,即总的taxa和在VCF文件中的taxa
                String[] taxainVCFfile = null;
                List<String> taxainVCFlist = new ArrayList();
                while ((temp = br.readLine()) != null) {

                    if (temp.startsWith("##")) {
                        continue;
                    }
                    if (temp.startsWith("#CHROM")) {//说明进入taxa信息列
                        List<String> l = PStringUtils.fastSplit(temp, "\t");
                        taxainVCFfile = new String[l.size() - 9];
                        for (int i = 9; i < l.size(); i++) {
                            taxainVCFfile[i - 9] = l.get(i);
                            hmtaxainVCFindex.put(taxainVCFfile[i - 9], i); //第0个taxa的genotype在第9行，第一个taxa的genotype在第10行，依次类推；
                        }

                        // ************** 在总的taxa中，搜索VCF文件中的taxa的index
                        for (int i = 0; i < taxainVCFfile.length; i++) { //找到taxa
                            int index = Arrays.binarySearch(taxa, taxainVCFfile[i]);
                            if (index > -1) {
                                hmtaxainTaxaindex.put(taxainVCFfile[i], index);
                            }
                        }
                        taxainVCFlist = Arrays.asList(taxainVCFfile);
                    }
                    if (!temp.startsWith("#")) {
                        cnt++;
                        if (cnt % 1000 == 0) {
                            System.out.println(String.valueOf(cnt) + " lines on chr " + String.valueOf(chr));
                        }

                        List<String> l = PStringUtils.fastSplit(temp.substring(0, 100), "\t");
                        int pos = Integer.valueOf(l.get(1));
                        int index = Arrays.binarySearch(delePos[chrIndex], pos);
                        if (index < 0) {
                            continue; //前面建立的 delePos[][] 和deleChar[][] 都是为现在在vcf文件中找位置贡献的，不是有害突变的位点，都过滤。*************************************************************
                        }
                        l = PStringUtils.fastSplit(temp, "\t");
                        int[] idx = new int[2];
                        if (l.get(3).charAt(0) == deleChar[chrIndex][index]) { //如果ref allele = deleChar allele
                            idx[0] = 0; //等于0则代表是 derived allele
                            idx[1] = 1;
                        } else {
                            idx[0] = 1;
                            idx[1] = 0;
                        }


                        //合计642个taxa，在A Bgenome中只有606（419+187）个，在Dgenome中只有455（419+36）个，我们要找到每个VCF文件中的genotype所对应的taxa的index
                        //code:
                        for (int i = 0; i < taxainVCFlist.size(); i++) {
                            int genotypeIndex = hmtaxainVCFindex.get(taxainVCFlist.get(i)); //获取该taxa的基因型所在的列的索引
                            int taxaIndex = hmtaxainTaxaindex.get(taxainVCFlist.get(i)); //获取该taxa所在总的642个数组中的索引，为后续写文件进行统计
                            String genoS = l.get(genotypeIndex); //指的是 GT:AD:GL 信息
                            if (genoS.startsWith(".")) {
                                continue; //如果以.开头，说明没有基因型信息，此位点没有测到。
                            }
                            List<String> ll = PStringUtils.fastSplit(genoS, ":"); //分开为GT   AD   GL三类

                            List<String> lll = PStringUtils.fastSplit(ll.get(1), ","); //lll指将AD提取出来，并以"，"号分割。如 0/0:1,0:0,3,28中 ，1，0分别代表ref和alt的测序深度
                            int depth = Integer.valueOf(lll.get(0)) + Integer.valueOf(lll.get(1)); //总得测序深度等于 ref + alt
                            if (depth < minDepth) {
                                continue; //最小测序深度是2，如果小于2，则弃用
                            }
                            lll = PStringUtils.fastSplit(ll.get(0), "/"); //这里lll指的是基因型GT，lll被重新赋值，之前代表的是AD
                            int v1 = Integer.valueOf(lll.get(0)); //v1等于 ref
                            int v2 = Integer.valueOf(lll.get(1)); // v2 等于 alt
                            int sum = 0;
                            //如果ref是derived allele，那么idx[0]=0;idx[1]=1. 当时0/0时，sum=2; 0/1时，sum=1; 1/1时， sum=0.
                            //如果alt是derived allele，那么idx[0]=1;idx[1]=0. 当时0/0时，sum=0; 0/1时，sum=1; 1/1时， sum=2.
                            //
                            //
                            if (v1 == idx[0]) { //
                                sum++;
                            }
                            if (v2 == idx[0]) {
                                sum++;
                            }
                            if (sum == 0) {
                            } else if (sum == 1) {
                                addCount[taxaIndex] += 0.5;
                            } else {
                                addCount[taxaIndex] += 1;
                                recCount[taxaIndex] += 1;
                            }
                            siteWithMinDepthCount[taxaIndex]++; //taxa有多少个有害突变位点
                        }

                    }

                }
                br.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });

        try {
            BufferedWriter bw = IOUtils.getTextWriter(addCountFileS);
            bw.write("Taxa\tDeleteriousCountPerHaplotype\tSiteCountWithMinDepth"); //每个taxa有多少个加性效应的derivedAllele 每个taxa在del库中含有基因型的个数
            bw.newLine();
            for (int i = 0; i < addCount.length; i++) {
                bw.write(taxa[i] + "\t" + String.valueOf(addCount[i]) + "\t" + String.valueOf(siteWithMinDepthCount[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
            bw = IOUtils.getTextWriter(recCountFileS);
            bw.write("Taxa\tDeleteriousCountPerLine\tSiteCountWithMinDepth");
            bw.newLine();
            for (int i = 0; i < recCount.length; i++) {
                bw.write(taxa[i] + "\t" + String.valueOf(recCount[i]) + "\t" + String.valueOf(siteWithMinDepthCount[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * 最开始的genic anno库进行的测试，没有添加 derivedSIFT列的信息
     */
    public void countDeleteriousVMapII_old() {
        String delVCFDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/002_vcf/del"; //有害变异的VCF文件路径
        String deleFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/001_snp/del_merged/chr_SNP_anno_del.txt.gz"; //有害变异信息库
        String addCountFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/additiveDeleterious_vmap2.txt"; //有害变异加性模型输出文件
        String recCountFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/recessiveDeleterious_vmap2.txt"; //有害变异隐形模型输出文件
        int minDepth = 2;//inclusive
        int chrNum = 42;
        ArrayList<Integer> chrList = new ArrayList();
        for (int i = 0; i < chrNum; i++) {
            chrList.add(i + 1);
        }


        RowTable<String> t = new RowTable(deleFileS);
        TIntArrayList[] posList = new TIntArrayList[chrNum];
        TCharArrayList[] charList = new TCharArrayList[chrNum];
        for (int i = 0; i < chrNum; i++) { //集合类数组，要初始化每一个list
            posList[i] = new TIntArrayList();
            charList[i] = new TCharArrayList();
        }
        new AoFile().readheader(deleFileS);
        String derivedAllele = null;
        for (int i = 0; i < t.getRowNumber(); i++) {
            int index = t.getCellAsInteger(i, 1) - 1; //染色体号的索引
            String ancestralAllele = t.getCell(i, 14);
            String majorAllele = t.getCell(i, 5);
            String minorAllele = t.getCell(i, 6);
            if (ancestralAllele.equals(majorAllele)) {
                derivedAllele = minorAllele;
                posList[index].add(t.getCellAsInteger(i, 2)); //将包含有derived allele的位点添加到Poslist
                charList[index].add(derivedAllele.charAt(0)); //返回的是char类型的字符
            }
            if (ancestralAllele.equals(minorAllele)) {
                derivedAllele = majorAllele;
                posList[index].add(t.getCellAsInteger(i, 2));
                charList[index].add(derivedAllele.charAt(0));
            }
        }
        int[][] delePos = new int[chrNum][];
        char[][] deleChar = new char[chrNum][];
        for (int i = 0; i < chrNum; i++) { //将每一个list转化为数组
            delePos[i] = posList[i].toArray();
            deleChar[i] = charList[i].toArray();
            Arrays.sort(delePos[i]);
        }


        String vmap2TaxaList = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/004_vmap2.1_TaxaList/TaxaList.txt";
        List<String> taxaList = new AoFile().getStringListwithoutHeader(vmap2TaxaList, 0);
        String[] taxa = taxaList.toArray(new String[taxaList.size()]);
        Arrays.sort(taxa);


        int taxaNum = taxa.length;
        double[] addCount = new double[taxa.length];
        int[] recCount = new int[taxa.length];
        int[] siteWithMinDepthCount = new int[taxa.length]; //每个taxa的深度

        chrList.parallelStream().forEach(chr -> {
            String delVmapFileS = "chr" + PStringUtils.getNDigitNumber(3, chr) + "_del_vmap2.1.vcf.gz";
            delVmapFileS = new File(delVCFDirS, delVmapFileS).getAbsolutePath();
            BufferedReader br = IOUtils.getTextGzipReader(delVmapFileS);
            int chrIndex = chr - 1;
            try {
                String temp = null;
                int cnt = 0;
                HashMap<String, Integer> hmtaxainVCFindex = new HashMap<>();
                HashMap<String, Integer> hmtaxainTaxaindex = new HashMap<>();
                //这里涉及2个数组概念， taxa and taxainVCFfile,即总的taxa和在VCF文件中的taxa
                String[] taxainVCFfile = null;
                List<String> taxainVCFlist = new ArrayList();
                while ((temp = br.readLine()) != null) {

                    if (temp.startsWith("##")) {
                        continue;
                    }
                    if (temp.startsWith("#CHROM")) {//说明进入taxa信息列
                        List<String> l = PStringUtils.fastSplit(temp, "\t");
                        taxainVCFfile = new String[l.size() - 9];
                        for (int i = 9; i < l.size(); i++) {
                            taxainVCFfile[i - 9] = l.get(i);
                            hmtaxainVCFindex.put(taxainVCFfile[i - 9], i); //第0个taxa的genotype在第9行，第一个taxa的genotype在第10行，依次类推；
                        }

                        // ************** 在总的taxa中，搜索VCF文件中的taxa的index
                        for (int i = 0; i < taxainVCFfile.length; i++) { //找到taxa
                            int index = Arrays.binarySearch(taxa, taxainVCFfile[i]);
                            if (index > -1) {
                                hmtaxainTaxaindex.put(taxainVCFfile[i], index);
                            }
                        }
                        taxainVCFlist = Arrays.asList(taxainVCFfile);
                    }
                    if (!temp.startsWith("#")) {
                        cnt++;
                        if (cnt % 1000 == 0) {
                            System.out.println(String.valueOf(cnt) + " lines on chr " + String.valueOf(chr));
                        }

                        List<String> l = PStringUtils.fastSplit(temp.substring(0, 100), "\t");
                        int pos = Integer.valueOf(l.get(1));
                        int index = Arrays.binarySearch(delePos[chrIndex], pos);
                        if (index < 0) {
                            continue; //前面建立的 delePos[][] 和deleChar[][] 都是为现在在vcf文件中找位置贡献的，不是有害突变的位点，都过滤。*************************************************************
                        }
                        l = PStringUtils.fastSplit(temp, "\t");
                        int[] idx = new int[2];
                        if (l.get(3).charAt(0) == deleChar[chrIndex][index]) { //如果ref allele = deleChar allele
                            idx[0] = 0; //等于0则代表是 derived allele
                            idx[1] = 1;
                        } else {
                            idx[0] = 1;
                            idx[1] = 0;
                        }


                        //合计642个taxa，在A Bgenome中只有606（419+187）个，在Dgenome中只有455（419+36）个，我们要找到每个VCF文件中的genotype所对应的taxa的index
                        //code:
                        for (int i = 0; i < taxainVCFlist.size(); i++) {
                            int genotypeIndex = hmtaxainVCFindex.get(taxainVCFlist.get(i)); //获取该taxa的基因型所在的列的索引
                            int taxaIndex = hmtaxainTaxaindex.get(taxainVCFlist.get(i)); //获取该taxa所在总的642个数组中的索引，为后续写文件进行统计
                            String genoS = l.get(genotypeIndex); //指的是 GT:AD:GL 信息
                            if (genoS.startsWith(".")) {
                                continue; //如果以.开头，说明没有基因型信息，此位点没有测到。
                            }
                            List<String> ll = PStringUtils.fastSplit(genoS, ":"); //分开为GT   AD   GL三类

                            List<String> lll = PStringUtils.fastSplit(ll.get(1), ","); //lll指将AD提取出来，并以"，"号分割。如 0/0:1,0:0,3,28中 ，1，0分别代表ref和alt的测序深度
                            int depth = Integer.valueOf(lll.get(0)) + Integer.valueOf(lll.get(1)); //总得测序深度等于 ref + alt
                            if (depth < minDepth) {
                                continue; //最小测序深度是2，如果小于2，则弃用
                            }
                            lll = PStringUtils.fastSplit(ll.get(0), "/"); //这里lll指的是基因型GT，lll被重新赋值，之前代表的是AD
                            int v1 = Integer.valueOf(lll.get(0)); //v1等于 ref
                            int v2 = Integer.valueOf(lll.get(1)); // v2 等于 alt
                            int sum = 0;
                            //如果ref是derived allele，那么idx[0]=0;idx[1]=1. 当时0/0时，sum=2; 0/1时，sum=1; 1/1时， sum=0.
                            //如果alt是derived allele，那么idx[0]=1;idx[1]=0. 当时0/0时，sum=0; 0/1时，sum=1; 1/1时， sum=2.
                            //
                            //
                            if (v1 == idx[0]) { //
                                sum++;
                            }
                            if (v2 == idx[0]) {
                                sum++;
                            }
                            if (sum == 0) {
                            } else if (sum == 1) {
                                addCount[taxaIndex] += 0.5;
                            } else {
                                addCount[taxaIndex] += 1;
                                recCount[taxaIndex] += 1;
                            }
                            siteWithMinDepthCount[taxaIndex]++; //taxa有多少个有害突变位点
                        }

                    }

                }
                br.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });

        try {
            BufferedWriter bw = IOUtils.getTextWriter(addCountFileS);
            bw.write("Taxa\tDeleteriousCountPerHaplotype\tSiteCountWithMinDepth"); //每个taxa有多少个加性效应的derivedAllele 每个taxa在del库中含有基因型的个数
            bw.newLine();
            for (int i = 0; i < addCount.length; i++) {
                bw.write(taxa[i] + "\t" + String.valueOf(addCount[i]) + "\t" + String.valueOf(siteWithMinDepthCount[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
            bw = IOUtils.getTextWriter(recCountFileS);
            bw.write("Taxa\tDeleteriousCountPerLine\tSiteCountWithMinDepth");
            bw.newLine();
            for (int i = 0; i < recCount.length; i++) {
                bw.write(taxa[i] + "\t" + String.valueOf(recCount[i]) + "\t" + String.valueOf(siteWithMinDepthCount[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void mergeDelgenicSNPAnnotation() {

        //change
//        String infileDirS = "";
//        String outfileS = "";
//        new CountSites().mergeTxtbysuffix(infileDirS, outfileS, ".txt.gz");


//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/001_snp/del";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/001_snp/del_merged/chr_SNP_anno_del.txt.gz";
//        new CountSites().mergeTxtbysuffix(infileDirS, outfileS, ".txt.gz");

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/001_snp/syn";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/001_snp/syn_merged/chr_SNP_anno_syn.txt.gz";
        new CountSites().mergeTxtbysuffix(infileDirS, outfileS, ".txt.gz");

        //Total lines without header count is 91306 at merged file /Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/001_snp/del_merged/chr_SNP_anno_del.txt.gz
        //12/26/19 8:37 PM
    }

}
