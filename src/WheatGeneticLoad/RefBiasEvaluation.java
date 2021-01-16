package WheatGeneticLoad;

import AoUtils.AoFile;
import AoUtils.CalVCF;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class RefBiasEvaluation {
    public RefBiasEvaluation(){
        this.extractAAF();
        this.checkAnnotationDBisinExonVCF();

    }

    public void GetSubspeciesDAF(){
        String infileS = "";

    }

    /**
     * 判断 annnotation DB 中的 site 是否都在 exon VCF 中
     * 结果： 全部都在！
     */

    public void checkAnnotationDBisinExonVCF(){

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/018_exonSNPAnnotation";
//        String aaffileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/016_exonVCF";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/001_ifAnnotationExonVCF";


        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/rscript/referenceEvaluation/data/002_popDepthBP_addDerivedSIFT";
        String aaffileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/rscript/referenceEvaluation/data/006_extractAAF_fromExonVCF";
        String outfileDirS = "/Users/Aoyue/Documents/test";

        List<File> fsList = AoFile.getFileListInDir(infileDirS);
        fsList.stream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                String aafFileS = new File(aaffileDirS,f.getName().replaceFirst("_SNP_anno.txt","_aaf.txt.gz")).getAbsolutePath();
                String outfileS = new File(outfileDirS, f.getName().split(".txt")[0] + ".txt").getAbsolutePath();

                // 现获取ExonVCF 的pos库
                TIntArrayList posList = AoFile.getTIntList(aafFileS,2);
                posList.sort();

                BufferedReader br = AoFile.readFile(infileS);
                BufferedWriter bw = AoFile.writeFile(outfileS);
                String header = br.readLine();
                String temp = null;
                List<String> l = new ArrayList<>();
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    String chr = l.get(1);
                    int pos = Integer.parseInt(l.get(2));
                    int index = posList.binarySearch(pos);
                    if (index < 0){
                        System.out.println(chr + "\t" + pos);
                    }
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


    public void extractAAF() {
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/002_exonSNPVCF";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/rscript/referenceEvaluation/data/006_extractAAF_fromExonVCF";

        List<File> fList = AoFile.getFileListInDir(infileDirS);
        Collections.sort(fList);

        String hexaFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/010_taxaList/BreadWheat_S420.txt";
        String diFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/010_taxaList/Ae.tauschii_S36.txt";
        String tetraFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/010_taxaList/EmmerWheat_S187.txt";

        String[] abdTaxa = AoFile.getStringArraybyList_withoutHeader(hexaFileS,0);
        String[] abTaxa = AoFile.getStringArraybyList_withoutHeader(tetraFileS,0);
        String[] dTaxa = AoFile.getStringArraybyList_withoutHeader(diFileS,0);
        Arrays.sort(abdTaxa);
        Arrays.sort(dTaxa);
        Arrays.sort(abTaxa);

        fList.parallelStream().forEach(f -> {
            String outfileS = new File(outfileDirS,f.getName().split("_exon")[0] + "_aaf.txt.gz").getAbsolutePath();
            int chr = Integer.parseInt(f.getName().split("_")[0].replaceFirst("chr", ""));
            String subgenome = RefV1Utils.getSubgenomeFromChrID(chr);


            List<Integer> indexABD = new ArrayList<>();
            List<Integer> indexABorD = new ArrayList<>();
            String[] taxaABorDArray = null;
            String aaf = null;  //在注释文件中是写AAF_D 还是AAF_AB

            if (subgenome.equals("D")) {
                taxaABorDArray = dTaxa;
                aaf = "AAF_D";
            } else {
                taxaABorDArray = abTaxa;
                aaf = "AAF_AB";
            }

            try {
                BufferedReader br = AoFile.readFile(f.getAbsolutePath());
                BufferedWriter bw = AoFile.writeFile(outfileS);
                bw.write("ID\tChr\tPos\tRef\tAlt\tMAF\tAAF_ABD\tAAF_AB");
                bw.newLine();
                String temp = null;
                List<String> l = new ArrayList<>();
                int cntSNP = 0; //totalSNP
                int cntkept = 0;

                while ((temp = br.readLine()) != null) {
                    //*********************** # section ************************************//
                    if (temp.startsWith("##")) continue;
                    if (temp.startsWith("#CHROM")) {
                        l = PStringUtils.fastSplit(temp);

                        for (int i = 9; i < l.size(); i++) {
                            String taxon = l.get(i);
                            int index1 = Arrays.binarySearch(abdTaxa, taxon);
                            int index2 = Arrays.binarySearch(taxaABorDArray, taxon);

                            if (index1 > -1) {
                                indexABD.add(i);
                            }
                            if (index2 > -1) {
                                indexABorD.add(i);
                            }
                        }
                        Collections.sort(indexABD);
                        Collections.sort(indexABorD);
                    }
                    //*********************** pos section ************************************//
                    if (!temp.startsWith("#")) { //
                        cntSNP++;
                        l = PStringUtils.fastSplit(temp);
                        String altList = l.get(4);
                        int pos  = Integer.parseInt(l.get(1));
                        String chrS = l.get(0);
                        String id = l.get(2);

                        List<String> lgeno = new ArrayList<>();
                        List<String> lABDGeno = new ArrayList<>();
                        List<String> lABorDGeno = new ArrayList<>();

                        for (int i = 9; i < l.size(); i++) {
                            lgeno.add(l.get(i));
                        }
                        for (int i = 0; i < indexABD.size(); i++) {
                            lABDGeno.add(l.get(indexABD.get(i)));
                        }
                        for (int i = 0; i < indexABorD.size(); i++) {
                            lABorDGeno.add(l.get(indexABorD.get(i)));
                        }

                        String[] genoArray = lgeno.toArray(new String[lgeno.size()]);
                        String[] hexaGenoArray = lABDGeno.toArray(new String[lABDGeno.size()]);
                        String[] ABorDGenoArray = lABorDGeno.toArray(new String[lABorDGeno.size()]);

                        String maf = this.getInfo(genoArray, altList);
                        String hexaAAF = this.getSubgenomeInfo(hexaGenoArray, altList).split(",")[0];
                        String ABorDAAF = this.getSubgenomeInfo(ABorDGenoArray, altList).split(",")[0];

                        StringBuilder sb = new StringBuilder();

                        sb.append(l.get(2)).append("\t").append(l.get(0)).append("\t").append(l.get(1)).append("\t").append(l.get(3)).append("\t").append(l.get(4)).append("\t");

                        sb.append(maf).append("\t").append(hexaAAF).append("\t").append(ABorDAAF);

                        bw.write(sb.toString());
                        bw.newLine();
                    }
                }
                br.close();
                bw.flush();
                bw.close();
                System.out.println( cntSNP + "\ttotal\t" + outfileS);

            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }

        });

    }


    private String getInfo(String[] genoArray, String altList) {
        int dp = 0; //总深度
        int nz = 0; //有基因型的个体数
        int nAlt = PStringUtils.fastSplit(altList, ",").size();
        int[] adCnt = new int[1 + nAlt]; //所有包括ref和alt的深度统计
        int[] acCnt = new int[1 + nAlt]; //所有包括ref和alt的基因型统计
        int[][] gnCnt = new int[1 + nAlt][1 + nAlt]; //GN到底代表什么？
        int ht = 0;
        List<String> tempList = null;
        List<String> temList = null;
        for (int i = 0; i < genoArray.length; i++) {
            if (genoArray[i].startsWith(".")) {
                nz++;
                continue;
            }
            tempList = PStringUtils.fastSplit(genoArray[i], ":"); //tempList是包含基因型AD还有PL的集合

            //先计算深度
            temList = PStringUtils.fastSplit(tempList.get(1), ","); //temList是AD所有的深度集合
            for (int j = 0; j < temList.size(); j++) {
                int c = Integer.parseInt(temList.get(j)); //c是第j个allele的深度值。注意AD的第一个是ref，第二个是次等位位点的深度，第三个是最小等位位点的深度
                dp += c; //dp是总深度
                adCnt[j] += c; //adCnt[j] 是第j个allele的深度值的总和，AD按照 ref alt1 alt2排序
            }

            //再计算基因型
            temList = PStringUtils.fastSplit(tempList.get(0), "/"); //temList是包含基因型拆分后的集合
            for (int j = 0; j < temList.size(); j++) { //0/0:13,0:0,4,25
                int c = Integer.parseInt(temList.get(j)); // c是基因型0 1 2 其中的一个
                acCnt[c]++; //acCnt[c] 是所有taxa基因型某一数值如0 1 2的总和
            }
            int index1 = Integer.parseInt(temList.get(0)); //0/0基因型的
            int index2 = Integer.parseInt(temList.get(1));
            gnCnt[index1][index2]++; //gnCnt[][]是二维数组，代表alt的个数的矩阵，比如有1个alt，则gnCnt[][]有gnCnt[0][0]  gnCnt[0][1] gnCnt[1][0] gnCnt[1][1]
            if (index1 != index2) {
                ht++;
            }
        }
        nz = genoArray.length - nz;
        int sum = 0; //所有allele的总数
        for (int i = 0; i < acCnt.length; i++) {
            sum += acCnt[i];
        }
        float maf = (float) ((double) acCnt[0] / sum);
        if (maf >= 0.5) {
            maf = (float) ((double) acCnt[1] / sum);
        }

        StringBuilder sb = new StringBuilder();

//        sb.append("DP=").append(dp).append(";NZ=").append(nz).append(";AD=");
//        for (int i = 0; i < adCnt.length; i++) {
//            sb.append(adCnt[i]).append(",");
//        }
//        sb.deleteCharAt(sb.length() - 1); //删除最后一个字符","号
//        sb.append(";AC=");
//        for (int i = 1; i < acCnt.length; i++) {
//            sb.append(acCnt[i]).append(",");
//        }
//        sb.deleteCharAt(sb.length() - 1);
//        sb.append(";GN=");
//        for (int i = 0; i < gnCnt.length; i++) { //二维数组的长度是第一维的长度，这里是2（只有1个alt） 或者3 (有2个alt)
//            for (int j = i + 1; j < gnCnt.length; j++) {
//                sb.append(gnCnt[i][j]).append(",");
//            }
//        }
//        sb.deleteCharAt(sb.length() - 1);
//        sb.append(";HT=").append(ht).append(";MAF=").append(String.format("%.4f", maf));
        sb.append(String.format("%.4f", maf));
        return sb.toString();
    }

    private String getSubgenomeInfo(String[] PopGenoArray, String altList) {
//        int   dp = 0;
        int nz = 0;
        int nAlt = PStringUtils.fastSplit(altList, ",").size();
//        int[] adCnt = new int[1 + nAlt]; //所有包括ref和alt的个数
        int[] acCnt = new int[1 + nAlt]; //所有包括ref和alt的个数
        int[][] gnCnt = new int[1 + nAlt][1 + nAlt]; //GN到底代表什么？
//        int ht = 0;
        List<String> tempList = null;
        List<String> temList = null;
        for (int i = 0; i < PopGenoArray.length; i++) {
            if (PopGenoArray[i].startsWith(".")) {
                nz++;
                continue;
            }
            tempList = PStringUtils.fastSplit(PopGenoArray[i], ":"); //tempList是包含基因型AD还有PL的集合
//            temList = PStringUtils.fastSplit(tempList.get(1), ","); //temList是AD所有的深度集合
//            for (int j = 0; j < temList.size(); j++) {
//                int c = Integer.parseInt(temList.get(j)); //c是第j个allele的深度值。注意AD的第一个是ref，第二个是次等位位点的深度，第三个是最小等位位点的深度
//                dp += c; //dp是总深度
//                adCnt[j] += c; //adCnt[j] 是第j个allele的深度值的总和，AD按照 ref alt1 alt2排序
//            }
            temList = PStringUtils.fastSplit(tempList.get(0), "/"); //temList是包含基因型拆分后的集合
            for (int j = 0; j < temList.size(); j++) {
                int c = Integer.parseInt(temList.get(j)); // c是基因型第j个数值
                acCnt[c]++; //acCnt[c] 是所有taxa基因型某一数值如0 1 2的总和
            }
//            int index1 = Integer.parseInt(temList.get(0)); //
//            int index2 = Integer.parseInt(temList.get(1));
//            gnCnt[index1][index2]++; //gnCnt[][]是二维数组，代表alt的个数的矩阵，比如有1个alt，则gnCnt[][]有gnCnt[0][0]  gnCnt[0][1] gnCnt[1][0] gnCnt[1][1]
//            if (index1 != index2) {
//                ht++;
//            }
        }
//        nz = PopGenoArray.length - nz;
        float missRate = (float) ((double) nz/PopGenoArray.length);

        int sum = 0;
        for (int i = 0; i < acCnt.length; i++) {
            sum += acCnt[i];
        }
        float maf = (float) ((double) acCnt[0] / sum);
        if (maf >= 0.5) {
            maf = (float) ((double) acCnt[1] / sum);
        }
        float aaf = (float) ((double) acCnt[1] / sum);

        StringBuilder sb = new StringBuilder();
        sb.append(String.format("%.4f", aaf)).append(",").append(String.format("%.4f", maf)).append(",").append(String.format("%.4f",missRate)); //.append(";MAF=")根据实际情况书写MAF_ABD MAF_AB MAF_D
        return sb.toString();
    }

}
