package WheatGeneticLoad;

import AoUtils.AoFile;
import AoUtils.AoMath;
import AoUtils.CalVCF;
import AoUtils.CountSites;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.table.RowTable;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;
import java.util.function.Predicate;

public class VMap2Cal {

    public VMap2Cal(){
//        this.extractAAFfromVMap2();
//        this.buildJointSFS();
        this.sampleSize2variantsDiscovery();

    }


    /**
     * 估算在exon区域，随着样本量的增大，各种类型的变异数目评估，判定是否达到饱和
     */
    public void sampleSize2variantsDiscovery(){
//        this.mergeExonVCF();
//        this.convert2GenoTable();
//        this.randomTaxa();
        this.getCount();


    }

    /**
     * 根据genotypetable，输出每种分类在每个亚基因组中的个数
     */
    public void getCount(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/021_exon_genoTable/test.txt";
//        List<String> l = Arrays.asList("2", "9", "0","2","2","9");
//        String test = this.genotypeCount(l);
//        System.out.println(test);



    }



    public boolean ifSegragation_ignoreIfRefisAnc(String genoClassCount){
        boolean out = false;
        List<String> l = PStringUtils.fastSplit(genoClassCount,";");
        int geno0 = Integer.parseInt(l.get(0).split("=")[1]);
        int geno1 = Integer.parseInt(l.get(1).split("=")[1]);
        int geno2 = Integer.parseInt(l.get(2).split("=")[1]);
        int geno9 = Integer.parseInt(l.get(3).split("=")[1]);
        if (geno1 >0 || geno2 >0){
            out = true;
        }else out = false;

        return out;
    }

    public boolean ifSegragationbyIfRefisAnc(String genoClassCount,String IfrefisAnc){
        boolean out = false;
        List<String> l = PStringUtils.fastSplit(genoClassCount,";");
        int geno0 = Integer.parseInt(l.get(0).split("=")[1]);
        int geno1 = Integer.parseInt(l.get(1).split("=")[1]);
        int geno2 = Integer.parseInt(l.get(2).split("=")[1]);
        int geno9 = Integer.parseInt(l.get(3).split("=")[1]);
        if (IfrefisAnc.equals("Anc")){
            if (geno1 >0 || geno2 >0){
                out = true;
            }else out = false;
        }
        if (IfrefisAnc.equals("Der")){
            if (geno0 >0 || geno1 >0){
                out = true;
            }else out = false;
        }
        return out;
    }

    /**
     * 判断一堆基因型中的 0/0 1/1 0/1 ./.个数
     * @return
     */
    public String genotypeCount (List<String> geno){
        String out = null;
        String[] genoClass = {"0","1","2","9"};
        int[] genoClassCount = new int[genoClass.length];
        Arrays.sort(genoClass);
        for (int i = 0; i < geno.size(); i++) {
            String query = geno.get(i);
            int index = Arrays.binarySearch(genoClass,query);
            if (index <0) System.out.println("There is some error about the genotype, please check your genotype.");
            genoClassCount[index]++;
        }

        StringBuilder sb = new StringBuilder();
        sb.append("0=").append(genoClassCount[0]).append(";").append("1=").append(genoClassCount[1]).append(";").
                append("2=").append(genoClassCount[2]).append(";").append("9=").append(genoClassCount[3]);
        out = sb.toString();
        return out;
    }

    public void randomTaxa(){
        //第一阶段： 抽样taxa，非连  续抽样
        String taxaListFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/011_taxaInfoDB/000_Dsub.txt";
        String[] taxaArray = AoFile.getStringArraybyList(taxaListFileS,0);
//        int[] sampleArray = {4,9};
        TIntArrayList sampleCountList = new TIntArrayList();
        sampleCountList.add(1);
        for (int i = 1; i < taxaArray.length; i++) { //每隔100个抽一次样******** 合计抽了6次 1，100，200，300，400，456
            int element = 100*i;
            if (element < taxaArray.length){
                sampleCountList.add(element);
            }else {
                sampleCountList.add(taxaArray.length);
                break;
            }
        }
        List<String>[] taxaListArray =  AoMath.continuousRandom(taxaArray,sampleCountList);


        //第二阶段：根据第一阶段抽样的taxa,进行genotype table 的提取，并返回数组类型的List
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/021_exon_genoTable/001_exon_Dsubgenome_genoTable.txt.gz";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/021_exon_genoTable/test.txt";
//        String[] taxaArray = {"BaiMaZha","BaiQiuMai"};
//        CalVCF.extractGenotable(infileS,taxaArray,outfileS);

        List<String> taxaList = new ArrayList<>(); taxaList.add("BaiMaZha");taxaList.add("BaiQiuMai");
        CalVCF.extractGenotable(infileS,taxaList,outfileS);


        //第三阶段：根据geontypeList,进行各种类型的计数





    }


    /**
     * 将exonVCF 转化为 derived allele table
     * 假设 ancestral allele 是 0/0 -> 0, derived allele 是1/1 -> 2, heter 是 0/1 ->2
     * coode: 0: ancestral allele 1:derived allele 2:derived
     */
    public void convert2GenoTable(){
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/020_exonVCF/ABsubgenome.vcf";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/021_exon_genoTable/001_exon_ABsubgenome_genoTable.txt.gz";
//        String taxaListFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/011_taxaInfoDB/000_ABsub.txt";
//        String[] taxaArray = AoFile.getStringArraybyList(taxaListFileS,0);
//        System.out.println(taxaArray.length);
//        CalVCF.extractVCFtable(infileS,taxaArray,outfileS);


        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/020_exonVCF/Dsubgenome.vcf";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/021_exon_genoTable/001_exon_Dsubgenome_genoTable.txt.gz";
        String taxaListFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/011_taxaInfoDB/000_Dsub.txt";
        String[] taxaArray = AoFile.getStringArraybyList(taxaListFileS,0);
        System.out.println(taxaArray.length);
        CalVCF.extractVCFtable(infileS,taxaArray,outfileS);

    }

    public void mergeExonVCF(){
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/016_exonVCF";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/020_exonVCF/001_exon.vcf";
//        CalVCF.mergeVCF(infileDirS,outfileS);


        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/016_exonVCF";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/020_exonVCF/";

        CountSites.mergeVCFtoAB_Dsubgenome(infileDirS,outfileDirS);

    }

    /**
     * 根据R 导入的数据，进行矩阵建立，最终画出2D-SFS图形
     */
    public void buildJointSFS(){

        int[][] matrix = new int[420*2][187*2];
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/rscript/referenceEvaluation/data/008_dafCount.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/rscript/referenceEvaluation/data/009_2dSFS.txt";
        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String header = br.readLine();
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            System.out.println(matrix.length);
            System.out.println(matrix[0].length);
            System.out.println();
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                int abd = Integer.parseInt(l.get(0));
                int ab = Integer.parseInt(l.get(1));
                matrix[abd-1][ab-1]++;
            }
            br.close();


            bw.write("\t");
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < matrix[0].length; i++) {
                sb.append("AB").append(i+1).append("\t");
            }
            sb.deleteCharAt(sb.length()-1);
            bw.write(sb.toString());
            bw.newLine();

            for (int i = 0; i < matrix.length; i++) {
                bw.write("ABD" + (i+1));
                for (int j = 0; j < matrix[0].length; j++) {
                    int count = matrix[i][j];
                    bw.write("\t" + count);
                }
                bw.newLine();
            }

            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    /**
     * Goal: 从VMap2.1 文件中提取第 8 列 （index 7），并判断maf值的大小
     */
    public void extractAAFfromVMap2(){
        //zcat chr1A_vmap2.1_hexaploid.vcf.gz |head -n 100 |grep -v "#" | awk '{print $8}'|cut -d";" -f7|cut -d"=" -f2| awk '{if ($0 <= 0.05) print $0}'|wc -l

//                String[] chrArr = {"1A", "1B", "2A", "2B", "3A", "3B", "4A", "4B", "5A", "5B", "6A", "6B", "7A", "7B"};
        String[] chrArr = {"1D","2D", "3D", "4D", "5D", "6D","7D"};
//        String[] chrArr = {"1A", "2A", "3A","4A", "5A", "6A", "7A", "1B", "2B", "3B", "4B", "5B", "6B", "7B", "1D", "2D", "3D",  "4D", "5D", "6D", "7D"};
//        String[] chrArr = {"1A", "2A", "3A","4A", "5A", "6A", "7A", "1B", "2B", "3B", "4B", "5B", "6B", "7B"};

//        String[] chrArr = {"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042"};
//        String[] chrArr ={"001","002","003","004","007","008","009","010","013","014","015","016","019","020","021","022","025","026","027","028","031","032","033","034","037","038","039","040"};
//        String[] chrArr ={"005","006","011","012","017","018","023","024","029","030","035","036","041","042"};

//        int[] darray = {5,6,11,12,17,18,23,24,29,30,35,36,41,42};

        for (int i = 0; i < chrArr.length; i++) {
            String chr = chrArr[i];
//            System.out.println( "zcat chr" + chr + "_vmap2.1_hexaploid.vcf.gz|grep -v \"#\" | awk '{print $8}'|cut -d\";\" -f7|cut -d\"=\" -f2| awk '{if ($0 <= 0.05) print $0}'|wc -l > ../test/" + chr + "_lessthan0.05.txt 2>&1 &");
//            System.out.println( "zcat chr" + chr + "_vmap2.1_hexaploid.vcf.gz|grep -v \"#\" | awk '{print $8}'|cut -d\";\" -f7|cut -d\"=\" -f2| awk '{if ($0 > 0.05) print $0}'|wc -l > ../test2/" + chr + "_morethan0.05.txt 2>&1 &");

//            System.out.println( "zcat chr" + chr + "_vmap2.1_tetraploid.vcf.gz|grep -v \"#\" | awk '{print $8}'|cut -d\";\" -f7|cut -d\"=\" -f2| awk '{if ($0 <= 0.05) print $0}'|wc -l > ../test3/" + chr + "_lessthan0.05.txt 2>&1 &");
//            System.out.println( "zcat chr" + chr + "_vmap2.1_tetraploid.vcf.gz|grep -v \"#\" | awk '{print $8}'|cut -d\";\" -f7|cut -d\"=\" -f2| awk '{if ($0 > 0.05) print $0}'|wc -l > ../test4/" + chr + "_morethan0.05.txt 2>&1 &");


//            System.out.println( "zcat chr" + chr + "_vmap2.1_diploid.vcf.gz|grep -v \"#\" | awk '{print $8}'|cut -d\";\" -f7|cut -d\"=\" -f2| awk '{if ($0 <= 0.05) print $0}'|wc -l > ../test5/" + chr + "_lessthan0.05.txt 2>&1 &");
            System.out.println( "zcat chr" + chr + "_vmap2.1_diploid.vcf.gz|grep -v \"#\" | awk '{print $8}'|cut -d\";\" -f7|cut -d\"=\" -f2| awk '{if ($0 > 0.05) print $0}'|wc -l > ../test6/" + chr + "_morethan0.05.txt 2>&1 &");

        }
    }

}
