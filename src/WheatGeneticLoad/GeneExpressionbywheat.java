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
        this.txtSinglethread();

    }

    /**
     * 选取CS的800个基因，进行相关性分析
     * 去除不表达的基因
     */
    public void txtSinglethread(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/001_geneExpressionDB/001_geneSummary_hexaploid_triad_Removed169_addTPM.txt.gz";
        String csgeneS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/002_merge/CS_geneSummary_triads_Remove169_wideTable.txt.gz";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/109_geneExpression/001_geneExpressionDB/002_geneSummary_CS_triad_Removed169_addTPM_removeNoexpression.txt.gz";

        AoFile.readheader(csgeneS);
        List<String> geneA = AoFile.getStringList(csgeneS,5);
        List<String> geneB = AoFile.getStringList(csgeneS,6);
        List<String> geneD = AoFile.getStringList(csgeneS,7);
        List<String> genesinCS = new ArrayList<>();
        genesinCS.addAll(geneA);genesinCS.addAll(geneB);genesinCS.addAll(geneD);

        try{
            BufferedReader br = AoFile.readFile(csgeneS);
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
                if (expressed.equals("FALSE"))continue;



            }
            br.close();

        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }


        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String header = br.readLine();
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                String gene = l.get(1).split("\\.")[0];

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
