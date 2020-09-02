package WheatGeneticLoad;

import AoUtils.AoFile;
import AoUtils.AoMath;
import AoUtils.AoString;
import pgl.infra.anno.gene.GeneFeature;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;

/**
 * @author AoyueBi
 * @data 2020-08-27 16:07
 */
public class AoWheatTriads {

    public AoWheatTriads(){
//        this.transformTriadsTable(); //for heatmap
//        this.getGeneInfo();
//        this.getDelHeter();
//        this.chromoMapInput();

    }


    public void chromoMapInput(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/004_eachGenewithPos/000_sourceData";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/004_eachGenewithPos/002";
        List<File> fsList = AoFile.getFileListInDir(infileDirS);
        AoFile.readheader("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/004_eachGenewithPos/000_sourceData/C1.triadpos.txt.gz");
        System.out.println("hhh");

        String geneFeatureFileS = "/Users/Aoyue/Documents/Data/wheat/gene/v1.1/wheat_v1.1_Lulab.pgf";
        GeneFeature gf = new GeneFeature(geneFeatureFileS);
        gf.sortGeneByName();
        fsList.parallelStream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                String outfileS = new File(outfileDirS, f.getName().split(".txt")[0] + ".txt").getAbsolutePath();
                BufferedReader br = AoFile.readFile(infileS);
                BufferedWriter bw = AoFile.writeFile(outfileS);
                String header = br.readLine();
                String temp = null;
                List<String> l = new ArrayList<>();
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    String gene = l.get(2);
                    String count = l.get(15);
                    if (Integer.parseInt(count)==0) continue;
                    int geneindex = gf.getGeneIndex(gene);
                    int start=gf.getGeneStart(geneindex);
                    int end = gf.getGeneEnd(geneindex);
                    int chrID=gf.getGeneChromosome(geneindex);
                    int refStart= RefV1Utils.getPosOnChromosome(chrID, start);
                    int refEnd=RefV1Utils.getPosOnChromosome(chrID, end);
                    String chr = AoString.getChrFromGene(gene);
                    bw.write(gene + "\t" + chr + "\t" + refStart + "\t" + refEnd + "\t" + count); bw.newLine();
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

    /**
     * 检查杂合子在 1 2 3 4 5 6 7 分别含有多少个基因
     */
    public void getDelHeter(){
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/004_eachGenewithPos/001_/003_C1_heter.txt";
//        String infileS ="/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/004_eachGenewithPos/001_/001_del_heter_annotation.txt";
//        AoMath.countCaseInGroup(infileS,9);
//        AoMath.countCaseInGroup(infileS,1);

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/004_eachGenewithPos/002";
        List<File> fsList = AoFile.getFileListInDir(infileDirS);
        for (int i = 0; i < fsList.size(); i++) {
            String infileS = fsList.get(i).getAbsolutePath();
            AoMath.countCaseInGroup(infileS,4);
            System.out.println(fsList.get(i).getName() + " done");
        }
    }

    /**
     * 为了将 del中是有害突变的杂合子在染色体上展示出来，特意提取heter信息，并整合做所需的文件
     * gene chrom start end count
     *
     * 测试用！！！！！！！！1
     *
     * // 注意pgf文件是1-42格式的
     */
    public void getGeneInfo(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/004_eachGenewithPos/001_/001_del_heter_annotation.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/004_eachGenewithPos/001_/002_del_heter_annotation.txt";

//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/004_eachGenewithPos/001_/003_C1_heter.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/004_eachGenewithPos/001_/004_C1_heter.txt";

        String geneFeatureFileS = "/Users/Aoyue/Documents/Data/wheat/gene/v1.1/wheat_v1.1_Lulab.pgf";
        GeneFeature gf = new GeneFeature(geneFeatureFileS);
        gf.sortGeneByName();
        RowTable<String> t = new RowTable<>(infileS);

        try{
            BufferedWriter bw = AoFile.writeFile(outfileS);
            for (int i = 0; i < t.getRowNumber(); i++) {
                String gene = t.getCell(i,0);
                String count = t.getCell(i,1);
//                String count = t.getCell(i,9);
                int geneindex = gf.getGeneIndex(gene);
                int start=gf.getGeneStart(geneindex);
                int end = gf.getGeneEnd(geneindex);
                int chrID=gf.getGeneChromosome(geneindex);
                int refStart= RefV1Utils.getPosOnChromosome(chrID, start);
                int refEnd=RefV1Utils.getPosOnChromosome(chrID, end);
                String chr = AoString.getChrFromGene(gene);
                bw.write(gene + "\t" + chr + "\t" + refStart + "\t" + refEnd + "\t" + count); bw.newLine();
            }
            bw.flush();
            bw.close();
        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * merge to one file and add chr
     */
    public void step5(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/001_sourceData/004_addChr";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/001_sourceData/005_mergetoOnefile/001_triadPosA.del.txt";
        File[] fs = AoFile.getFileArrayInDir(infileDirS);
        Arrays.sort(fs);
        try{
            BufferedReader[] br = new BufferedReader[fs.length];
            BufferedWriter bw = AoFile.writeFile(outfileS);
            for (int i = 0; i < fs.length; i++) {
                File f = fs[i];
                String infileS = f.getAbsolutePath();
                br[i] = AoFile.readFile(infileS);
            }

            // 开始读文件，每个文件的第一个文件全读，第二个第三个开始只从第3列开始读。
            String[] temp = new String[fs.length];

            while ((temp[0] = br[0].readLine()) != null) {
                bw.write(temp[0]); //读完第一行
                for (int i = 1; i < fs.length; i++) {
                    temp[i] = br[i].readLine();
                    List<String> l = PStringUtils.fastSplit(temp[i]);
                    for (int j = 2; j < l.size() ; j++) { //意味着从第3列开始读
                        bw.write("\t" + l.get(j));
                    }
                }
                bw.newLine();
            }
            for (int i = 0; i < br.length; i++) {
                br[i].close();
            }
            bw.flush();bw.close();
        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 将每个文件添加一行染色体号信息
     */
    public void step4(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/001_sourceData/003_mergeByrowByChr";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/001_sourceData/004_addChr";
        File[] fs = AoFile.getFileArrayInDir(infileDirS);
        Arrays.sort(fs);
        try{
            for (int i = 0; i < fs.length; i++) {
                File f = fs[i];
                String chr = f.getName().substring(0,1);
                String infileS = f.getAbsolutePath();
                String outfileS = new File(outfileDirS,f.getName()).getAbsolutePath();
                BufferedReader br = AoFile.readFile(infileS);
                BufferedWriter bw = AoFile.writeFile(outfileS);

                String header = br.readLine();
                List<String> l = PStringUtils.fastSplit(header);
                bw.write("Chr\tChr2");
                for (int j = 2; j < l.size(); j++) {
                    bw.write("\t" + chr);
                }
                bw.newLine();
                bw.write(header);bw.newLine();

                String temp = null;
                while ((temp = br.readLine()) != null) {
                    bw.write(temp);bw.newLine();
                }


                bw.flush();bw.close();br.close();
            }

        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 所有以1开头 以2开头的文件，横向合并起来
     */
    public void step3(String infileDirS, String outfileDirS){
//        String infileDirS = "";
//        String outfileDirS = "";

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/001_sourceData/002_transform";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/001_sourceData/003_mergeByrowByChr";


//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/003_dataProcessingfrom002/003";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/003_dataProcessingfrom002/004";

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/003_dataProcessingfrom002/VrnA1/003";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/003_dataProcessingfrom002/VrnA1/004";


//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/003_dataProcessingfrom002/Q/003";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/003_dataProcessingfrom002/Q/004";


//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/003_dataProcessingfrom002/TtBtr/003";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/003_dataProcessingfrom002/TtBtr/004";

        File[] fsall = AoFile.getFileArrayInDir(infileDirS);

        Set<String> chrSet = new HashSet<>();
        for (int i = 0; i < fsall.length; i++) {
            File f = fsall[i];
            String c = f.getName().substring(0,1); //获取染色体号
            chrSet.add(c);
        }
        //获取含有的染色体号
        String[] chrArray = chrSet.toArray(new String[chrSet.size()]);
        Arrays.sort(chrArray);
        for (int i = 0; i < chrArray.length; i++) {
            String chr = chrArray[i];
            String outfileS = new File(outfileDirS,chr + ".txt").getAbsolutePath();
            File[] fs = IOUtils.listFilesStartsWith(fsall,chr);
            AoFile.mergeTxt_byFileArray(fs,outfileS);
        }
    }

    public void step2(String infileDirS, String outfileDirS){
//        String infileDirS = "";
//        String outfileDirS = "";

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/001_sourceData/001_split";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/001_sourceData/002_transform";

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/003_dataProcessingfrom002/002";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/003_dataProcessingfrom002/003";

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/003_dataProcessingfrom002/VrnA1/002";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/003_dataProcessingfrom002/VrnA1/003";

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/003_dataProcessingfrom002/Q/002";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/003_dataProcessingfrom002/Q/003";

//                String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/003_dataProcessingfrom002/TtBtr/002";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/003_dataProcessingfrom002/TtBtr/003";


        File[] fs = AoFile.getFileArrayInDir(infileDirS);
        Arrays.sort(fs);

        try {
            for (int i = 0; i < fs.length; i++) {
                String infileS = fs[i].getAbsolutePath();
                String outfileS = new File(outfileDirS,fs[i].getName()).getAbsolutePath();
                BufferedReader br = AoFile.readFile(infileS);
                BufferedWriter bw = AoFile.writeFile(outfileS);
                String sub = fs[i].getName().substring(1,2);
                System.out.println("dealing with " + sub);

                int column = AoFile.countFileColumnNumber(infileS); //获取列数
                HashMap<Integer,String>[] hmArray = new HashMap[column-3]; //根据列数，建立pos taxa 对应的 HashMap
                for (int j = 0; j < hmArray.length; j++) {
                    hmArray[j] = new HashMap<>();
                }
                String[] taxaArray = new String[hmArray.length]; //建立taxa数组
                List<String> taxaList = new ArrayList<>(); //建立taxa集合
                List<String> headerList = AoFile.getHeader(infileS); //获取header集合
                for (int j = 3; j < headerList.size(); j++) {
                    taxaList.add(headerList.get(j).replaceFirst(".MeanNormalizedDeleteriousLoad",""));
                }
                taxaArray = taxaList.toArray(new String[taxaList.size()]);

                String temp = null;
                String header = br.readLine();
                List<String> l = new ArrayList<>();
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    cnt++;
                    int pos = Integer.parseInt(l.get(1));
                    for (int j = 3; j < l.size(); j++) {
                        String value = l.get(j);
                        hmArray[j-3].put(pos,value);
                    }
                }
                br.close();


                //开始写文件
                //写表头
                List<Integer> lpos = new ArrayList<>(hmArray[0].keySet());
//                Collections.sort(lpos);
                bw.write("Sub\tTaxa");
                for (int j = 0; j < lpos.size(); j++) {
                    bw.write("\t" + lpos.get(j));
                }
                bw.newLine();

                //写正文
                for (int j = 0; j < taxaArray.length; j++) {
                    String taxa = taxaArray[j];
                    bw.write(sub + "\t" + taxa);
                    for (int k = 0; k < lpos.size(); k++) { //数据转置
                        String test = hmArray[j].get(k);
                        bw.write("\t" + hmArray[j].get(lpos.get(k)));
                    }
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                System.out.println();
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }


    public void step1(String infileS, String outfileDirS){

//        String infileS = "";
//        String outfileDirS = "";

        try {
//            String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/001_sourceData/001_triadPosA.IndividualMerged.txt";
//            String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/001_sourceData/001_split";

//            String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/003_dataProcessingfrom002/001/PpdA1.txt";
//            String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/003_dataProcessingfrom002/002";

//            String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/003_dataProcessingfrom002/VrnA1/001/VrnA1.txt";
//            String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/003_dataProcessingfrom002/VrnA1/002";

//                    String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/003_dataProcessingfrom002/Q/001/Q.txt";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/003_dataProcessingfrom002/Q/002";

//                    String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/003_dataProcessingfrom002/TtBtr/001/TtBtrA1A2.txt";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/003_dataProcessingfrom002/TtBtr/002";

            List<String> lchr = AoFile.getStringListbySet(infileS,0);
            String[] chrArray = lchr.toArray(new String[lchr.size()]);
            Arrays.sort(chrArray);
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter[] bws = new BufferedWriter[chrArray.length];
            for (int i = 0; i < chrArray.length; i++) {
                String chr = chrArray[i];
                String outfileS = new File(outfileDirS,chr + ".txt").getAbsolutePath();
                bws[i] = AoFile.writeFile(outfileS);
            }
            String temp = null;
            String header = br.readLine();
            for (int i = 0; i < bws.length; i++) {
                bws[i].write(header); bws[i].newLine();
            }
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                String chr = l.get(0);
                int index = Arrays.binarySearch(chrArray,chr);
                bws[index].write(temp);bws[index].newLine();
            }
            br.close();
            for (int i = 0; i < bws.length; i++) {
                bws[i].flush();
                bws[i].close();
            }

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public void transformTriadsTable(){

//        String infileS = "";
//        String outfileDirS = "";

        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/002_ImportantGeneBlock/Q.txt";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/037_Triads/003_dataProcessingfrom002/test";
        String[] outfileSS = new String[3];
        for (int i = 0; i < outfileSS.length; i++) {
            String dir = PStringUtils.getNDigitNumber(3,i+1);
            outfileSS[i] = new File(outfileDirS,dir).getAbsolutePath();
            new File(outfileSS[i]).mkdirs();
        }

        this.step1(infileS,outfileSS[0]); //将文件按照染色体一条一条拆分
        this.step2(outfileSS[0],outfileSS[1]); //转置，并只取起始位置
        this.step3(outfileSS[1],outfileSS[2]); //将文件根据染色体号进行合并，自动判断有几条染色体就合并成几个文件

//        this.step4();
//        this.step5();
    }

}
