package WheatGeneticLoad;

import AoUtils.AoFile;
import AoUtils.Triads.Triadsgenes;
import org.apache.commons.lang.ArrayUtils;
import pgl.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

public class GeneExpressionbywheat {
    public GeneExpressionbywheat(){
//        this.addTPM();
        this.getCSgeneSummary();

    }

    /**
     * 选取CS的800个基因，进行相关性分析
     * 去除不表达的基因,并添加 ternary分组和是否共线等信息
     */
    public void getCSgeneSummary(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/001_geneExpressionDB/001_geneSummary_hexaploid_triad_Removed169_addTPM.txt.gz";
        String csgeneFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/002_merge/CS_geneSummary_triads_Remove169_wideTable.txt.gz";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/001_geneExpressionDB/002_geneSummary_CS_triad_Removed169_addTPM_removeNoexpression_removeNOHGDel.txt.gz";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/001_geneExpressionDB/002_geneSummary_CS_triad_Removed169_addTPM_removeNoexpression.txt.gz";

        List<String> genesinCS = new ArrayList<>();
        HashMap<String,String> hmTriadIDLoadGroup = new HashMap<>();
        HashMap<String,String> hmgeneIfexpressed = new HashMap<>();
        HashMap<String,String> hmgeneIfsyntenic = new HashMap<>();

        AoFile.readheader(csgeneFileS);
        List<String> geneA = AoFile.getStringList(csgeneFileS,5);
        List<String> geneB = AoFile.getStringList(csgeneFileS,6);
        List<String> geneD = AoFile.getStringList(csgeneFileS,7);

//        genesinCS.addAll(geneA);genesinCS.addAll(geneB);genesinCS.addAll(geneD);
//        HashMap<String,String> hmTriadIDLoadGroup = AoFile.getHashMapStringKey(csgeneFileS,0,4);

        try{
            BufferedReader br = AoFile.readFile(csgeneFileS);
            String header = br.readLine();
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                String triadID = l.get(0);
                String loadGroup = l.get(4);
                String genea = l.get(5);
                String geneb = l.get(6);
                String gened = l.get(7);
                String synteny = l.get(8);
                String expressed = l.get(9);

                /**
                 * 过滤NA 即没有 high del
                 */
                String contribution = l.get(1);
                if (contribution.equals("NA")) continue;

//                if (expressed.equals("FALSE"))continue;
                genesinCS.add(genea);genesinCS.add(geneb);genesinCS.add(gened);
                hmgeneIfexpressed.put(genea,expressed);
                hmgeneIfexpressed.put(geneb,expressed);
                hmgeneIfexpressed.put(gened,expressed);
                hmTriadIDLoadGroup.put(genea,loadGroup);
                hmTriadIDLoadGroup.put(geneb,loadGroup);
                hmTriadIDLoadGroup.put(gened,loadGroup);
                hmgeneIfsyntenic.put(genea,synteny);
                hmgeneIfsyntenic.put(geneb,synteny);
                hmgeneIfsyntenic.put(gened,synteny);
            }
            br.close();

        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        Collections.sort(genesinCS);
        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String header = br.readLine();
            bw.write(header + "\tTernaryGroup\tIfSyntenic");
            bw.newLine();
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            String syntenicState = null;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                String gene = l.get(1).split("\\.")[0];
                int index = Collections.binarySearch(genesinCS,gene);
                if (index < 0)continue;
                /**
                 * new method
                 */
                String express = hmgeneIfexpressed.get(gene);
                if (express.equals("FALSE"))continue;
                String loadGroup = hmTriadIDLoadGroup.get(gene);
                syntenicState = hmgeneIfsyntenic.get(gene);

//                boolean out = Triadsgenes.ifExpressedBasedGene(gene); //判断是否表达
//                if (!out)continue;
//                String triadID = Triadsgenes.getTriadID(gene);
//                String loadGroup = hmTriadIDLoadGroup.get(triadID); //获取 分组信息
//                boolean syntenic = Triadsgenes.ifSyntenicBasedTriadID(triadID);
//                if (syntenic) syntenicState = "1";
//                else syntenicState = "0";
                bw.write(temp+"\t"+ loadGroup + "\t" + syntenicState);
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

    /**
     * 向基因summary文件中添加 TPM 信息，即每个基因的平均表达量，最大值，最小值，标准误
     *
     * step1:
     */
    public  void addTPM (){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/002_merge/002_geneSummary_hexaploid_triad_Removed169.txt.gz";
        String geneExpressionS = "/Users/Aoyue/Documents/Data/wheat/gene/gene_expression_dontMove/geneExpression.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/001_geneExpressionDB/001_geneSummary_hexaploid_triad_Removed169_addTPM.txt.gz";

        AoFile.readheader(infileS);
        //建库HashMap
        HashMap<String,Double> hmgeneMean = AoFile.getHashMapdoubleValue(geneExpressionS,0,5);
        HashMap<String,Double> hmgeneSd = AoFile.getHashMapdoubleValue(geneExpressionS,0,6);
        HashMap<String,Double> hmgeneRsd = AoFile.getHashMapdoubleValue(geneExpressionS,0,7);

        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String header = br.readLine();
            bw.write(header + "\tTPM_mean\tTPM_sd\tTPM_rsd");
            bw.newLine();
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                String gene = l.get(1).split("\\.")[0];
                double meanTPM = hmgeneMean.get(gene);
                double sdTPM = hmgeneSd.get(gene);
                double rsdTPM = hmgeneRsd.get(gene);
                bw.write(temp+"\t"+String.format("%.4f",meanTPM) + "\t" + String.format("%.4f",sdTPM)+ "\t" + String.format("%.4f",rsdTPM));
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


}
