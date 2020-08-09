package PopulationAnalysis;


import AoUtils.*;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;

/**
 *
 * @author Aoyue
 */
public class AoHeterozygosity {
    public AoHeterozygosity(){
//        this.scriptSNPbased();
//        this.windowCal();
//        this.scriptforIndi();
//        this.script_calWindowStep();
//        this.mergeTxt();
//
//        /**
//         * 把个体在每个位点的基因型提取出来，并进行滑窗计算区段的杂合度，即在染色体2M内杂合子的个数除以2M的基因型个数
//         */
//        this.mkGenotype("/Users/Aoyue/Documents/ok/chr7B_vmap2.1_heter_SNPbased_Cultivar.vcf.gz","/Users/Aoyue/Documents/out/chr7B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz");
//        this.script_mkGenotype();
//
//        this.calWindowStep_RH_indivi(); //也是多线程
//        this.runJarParallele();
//        this.mergeTxt_calWindowStep_RH_indivi();

        /**
         * 新一波的计算
         */
//        this.getGenotypeTable();
//        this.changeChrPosandMergeGenoTable();
//        this.mkBinRH();
        this.mergeTxtandAddGroup();





    }


    public void mergeTxtandAddGroup() {
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/028_hapScannerAgain/008_indiviHeteralongChrom/004_mkBinTable";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/028_hapScannerAgain/008_indiviHeteralongChrom/005_mergefrom004/chr1A1B_Wild_emmer_RH_2000000_1000000.txt";
        File[] fs = AoFile.getFileArrayInDir(infileDirS);
        Arrays.sort(fs);
        try {
            String infileS = fs[0].getAbsolutePath();
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            /**
             * 需要改动,header的名字可以自定义
             */
            //read header
            bw.write(br.readLine() + "\tTaxa");
            bw.newLine();

            int cnttotal = 0;
            //read context
            for (int i = 0; i < fs.length; i++) {
                infileS = fs[i].getAbsolutePath();
                /**
                 * 需要改动
                 */
                String name = fs[i].getName();
                String group = name.split("emmer")[1].split("_RH")[0];

                br = AoFile.readFile(infileS);
                br.readLine(); //read header
                String temp = null;
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    cnttotal++;
                    StringBuilder sb = new StringBuilder();
                    sb.append(temp).append("\t").append(group);
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

    public void mkBinRH(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/028_hapScannerAgain/008_indiviHeteralongChrom/003_merge";
        String outfileDirS ="/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/028_hapScannerAgain/008_indiviHeteralongChrom/004_mkBinTable";
        int window = 2000000;
        int step = 1000000;

        String[] chrArr = {"1A", "1B"};
        for (int i = 0; i < chrArr.length; i++) {
            String infileS = new File(infileDirS,"chr" + chrArr[i] + "_Wild_emmer_genoTable.txt.gz").getAbsolutePath();
            List<String> headerList = AoFile.getheader(infileS);
            //CHROM	BIN_START	BIN_END	N_VARIANTS	HETEROZYGOSITY
            //CHROM	POS	PI362699	PI487260	PI94648	W11	W13	W14	W15	W16	W18	W19	W20	W5	W7	W8
            for (int j = 2; j < headerList.size(); j++) {
                String taxa = headerList.get(j);
                String outfileS = new File(outfileDirS,"chr" + chrArr[i] + "_Wild_emmer" + taxa + "_RH_" + window + "_" + step + ".txt").getAbsolutePath();
                Bin.ResidualHeterozygosity(infileS,chrArr[i],1,j,window,step,outfileS);
            }
            System.out.println();
        }
    }

    public void changeChrPosandMergeGenoTable(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/028_hapScannerAgain/008_indiviHeteralongChrom/002_out_genoTable";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/028_hapScannerAgain/008_indiviHeteralongChrom/003_merge";

        new CountSites().mergefileandChangeChrPos_chr1and2(infileDirS,outfileDirS);

    }



    /**
     * 获取相应taxa的基因型信息。结果为TXT格式的table
     */
    public void getGenotypeTable(){

        // just test
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/001_depth/001_subsetData/chr1D_vmap2_subset0.001.vcf.gz";
//        String taxaListFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/028_hapScannerAgain/008_indiviHeteralongChrom/001_test/Ae.tauschii.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/028_hapScannerAgain/008_indiviHeteralongChrom/001_test/chr1D_Ae.tauschii_genoTable.txt";
//        CalVCF.extractVCFtable(infileS,taxaListFileS,outfileS);

        // run HPC
        String infileDirS = "/data4/home/aoyue/vmap2/analysis/023_hapScanner_basedPopDepth/ab/out/VCF";
        String outfileDirS ="/data4/home/aoyue/vmap2/analysis/023_hapScanner_basedPopDepth/002_heter_Indiv/002_out_genoTable";
        String logDirS = "/data4/home/aoyue/vmap2/analysis/023_hapScanner_basedPopDepth/002_heter_Indiv/z_log";
        String taxaListFileS = "/data4/home/aoyue/vmap2/analysis/023_hapScanner_basedPopDepth/002_heter_Indiv/001_taxalist/Wild_emmer.txt";
        String name = new File(taxaListFileS).getName().split(".txt")[0];

        String[] chrArr ={"001","002","003","004"};
        for (int i = 0; i < chrArr.length; i++) {
            String infileS = new File(infileDirS,"chr" + chrArr[i] + ".vcf.gz").getAbsolutePath();
            String outfileS = new File(outfileDirS,"chr" + chrArr[i] + "_" + name + "_genoTable.txt.gz").getAbsolutePath();
            String logfileS = new File(logDirS,"log_" + new File(outfileS).getName().split(".gz")[0]).getAbsolutePath(); //不管是不是gz结尾，我们只取gz前的部分，妙！
            System.out.println("nohup java -jar 042_extractVCFtable.jar " + infileS + " " + taxaListFileS + " " + outfileS + " > " + logfileS  + " 2>&1 &" );
        }
    }


    public void mergeTxt_calWindowStep_RH_indivi(){
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/004_heterogozysity/004_indi_RH/001_out";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/004_heterogozysity/004_indi_RH/002_merge/Heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt";
//        new CountSites().mergeTxt(infileDirS,outfileS);

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/004_heterogozysity/004_indi_RH/003_out_DD";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/004_heterogozysity/004_indi_RH/004_merge/heter_SNPbased_Ae.tauschii_AE1211_chrposGenotype_RH_2Mwindow_1Mstep.txt";
//        new CountSites().mergeTxt(infileDirS,outfileS);

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/004_heterogozysity/004_indi_RH/006_out_AABB";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/004_heterogozysity/004_indi_RH/007_merge/heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype_RH_2Mwindow_1Mstep.txt";
        new CountSites().mergeTxt(infileDirS,outfileS);
    }


    public void runJarParallele(){
        String infileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans";
        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/DD";
        String logDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004";
//        String[] chrArr = {"1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B", "7D"};
        String[] chrArr = {"1A", "1B", "2A", "2B", "3A", "3B", "4A", "4B", "5A", "5B", "6A", "6B", "7A", "7B"};
//        String[] chrArr = {"1D","2D", "3D", "4D", "5D", "6D","7D"};
        for (int j = 0; j < chrArr.length; j++) {
//            String infileS = new File(infileDirS,"chr" + chrArr[j] + "_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz").getAbsolutePath();
//            String outfileS = new File(outfileDirS,"chr" + chrArr[j] +"_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt").getAbsolutePath();
//            String logfileS = new File(logDirS,"log_chr" + chrArr[j] + "_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt").getAbsolutePath();

//            String infileS = new File(infileDirS,"chr" + chrArr[j] + "_vmap2.1_heter_SNPbased_Ae.tauschii_AE1211_chrposGenotype.txt.gz").getAbsolutePath();
//            String outfileS = new File(outfileDirS,"chr" + chrArr[j] +"_vmap2.1_heter_SNPbased_Ae.tauschii_AE1211_chrposGenotype_RH_2Mwindow_1Mstep.txt").getAbsolutePath();
//            String logfileS = new File(logDirS,"log_chr" + chrArr[j] + "_vmap2.1_heter_SNPbased_Ae.tauschii_AE1211_chrposGenotype_RH_2Mwindow_1Mstep.txt").getAbsolutePath();

            String infileS = new File(infileDirS,"chr" + chrArr[j] + "_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype.txt.gz").getAbsolutePath();
            String outfileS = new File(outfileDirS,"chr" + chrArr[j] +"_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype_RH_2Mwindow_1Mstep.txt").getAbsolutePath();
            String logfileS = new File(logDirS,"log_" + new File(outfileS).getName().split(".gz")[0]).getAbsolutePath();


            System.out.println("nohup java -jar 035_calWindowStep_RH_indivi.jar " + infileS + " 2000000 1000000 " + outfileS + " > " + logfileS  + " 2>&1 &" );
        }

//        nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr1A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/chr1A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr1A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &
//        nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr2A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/chr2A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr2A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &
//        nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr2B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/chr2B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr2B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &
//        nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr3A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/chr3A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr3A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &
//        nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr3B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/chr3B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr3B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &
//        nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr3D_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/chr3D_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr3D_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &
//        nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr4A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/chr4A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr4A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &
//        nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr4B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/chr4B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr4B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &
//        nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr5A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/chr5A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr5A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &
//        nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr5B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/chr5B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr5B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &
//        nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr5D_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/chr5D_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr5D_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &
//        nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr6A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/chr6A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr6A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &
//        nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr6B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/chr6B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr6B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &
//        nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr7A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/chr7A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr7A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &
//        nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr7B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/chr7B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr7B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &
        //
        //nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr1D_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/chr1D_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr1D_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &
        //nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr2D_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/chr2D_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr2D_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &
        //nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr4D_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/chr4D_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr4D_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &
        //nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr6D_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/chr6D_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr6D_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &
        //nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr7D_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/chr7D_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr7D_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &

        //        nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr1B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/chr1B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr1B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &
//        nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr2A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/chr2A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr2A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &
//                nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr2B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/chr2B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr2B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &
//                nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr3B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/chr3B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr3B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &
//                nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr6B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/chr6B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr6B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &
//                nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr7A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/chr7A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr7A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &
//                nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr7B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/chr7B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr7B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &


        //nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr1A_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/DD/chr1A_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr1A_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &
        //nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr1B_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/DD/chr1B_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr1B_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &
        //nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr2A_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/DD/chr2A_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr2A_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &
        //nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr2B_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/DD/chr2B_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr2B_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &
        //nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr3A_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/DD/chr3A_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr3A_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &
        //nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr3B_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/DD/chr3B_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr3B_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &
        //nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr4A_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/DD/chr4A_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr4A_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &
        //nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr4B_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/DD/chr4B_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr4B_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &
        //nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr5A_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/DD/chr5A_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr5A_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &
        //nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr5B_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/DD/chr5B_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr5B_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &
        //nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr6A_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/DD/chr6A_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr6A_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &
        //nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr6B_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/DD/chr6B_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr6B_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &
        //nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr7A_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/DD/chr7A_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr7A_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &
        //nohup java -jar 035_calWindowStep_RH_indivi.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr7B_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype.txt.gz 2000000 1000000 /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/006_005_bin/DD/chr7B_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype_RH_2Mwindow_1Mstep.txt > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/004/log_chr7B_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype_RH_2Mwindow_1Mstep.txt 2>&1 &
    }



    /**
     * 根据点的数值，计算 滑窗的数值
     *
     *
     */
    public void calWindowStep_RH_indivi(){
        String infileDirS = "/Users/Aoyue/Documents/out";
        String outfileDirS = "/Users/Aoyue/Documents/out1";
        List<File> fsList = IOUtils.getVisibleFileListInDir(infileDirS);
        fsList.stream().forEach(f -> {
            String infileS = f.getAbsolutePath();
            String outfileS = null;
            if (infileS.endsWith(".txt")) {
                outfileS = new File(outfileDirS, f.getName().split(".txt")[0] + "_RH_2Mwindow_1Mstep.txt.gz").getAbsolutePath();
            } else if (infileS.endsWith(".txt.gz")) {
                outfileS = new File(outfileDirS, f.getName().split(".txt.gz")[0] + "_RH_2Mwindow_1Mstep.txt.gz").getAbsolutePath();
            }
            new Bin().calwindowstep_ResidualHeterozygosity(infileS,2000000,1000000,outfileS);
            System.out.println(f.getName() + "\tis completed at " + outfileS);
        });
    }


    public void script_mkGenotype(){
        String infileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/003_indiVCF";
        String outfileDirS ="/data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans";
        String logDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/003/";
//        String[] chrArr = {"1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B", "7D"};
        String[] chrArr = {"1A", "1B", "2A", "2B", "3A", "3B", "4A", "4B", "5A", "5B", "6A", "6B", "7A", "7B"};
//        String[] chrArr = {"1D","2D", "3D", "4D", "5D", "6D","7D"};
        for (int j = 0; j < chrArr.length; j++) {
//            String infileS = new File(infileDirS,"chr" + chrArr[j] + "_vmap2.1_heter_SNPbased_Cultivar.vcf.gz").getAbsolutePath();
//            String outfileS = new File(outfileDirS,"chr" + chrArr[j] + "_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz").getAbsolutePath();
//            String logfileS = new File(logDirS,"log_chr" + chrArr[j] + "_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt").getAbsolutePath();

//            String infileS = new File(infileDirS,"chr" + chrArr[j] + "_vmap2.1_heter_SNPbased_Ae.tauschii_AE1211.vcf.gz").getAbsolutePath();
//            String outfileS = new File(outfileDirS,"chr" + chrArr[j] + "_vmap2.1_heter_SNPbased_Ae.tauschii_AE1211_chrposGenotype.txt.gz").getAbsolutePath();
//            String logfileS = new File(logDirS,"log_chr" + chrArr[j] + "_vmap2.1_heter_SNPbased_Ae.tauschii_AE1211_chrposGenotype.txt").getAbsolutePath();

            String infileS = new File(infileDirS,"chr" + chrArr[j] + "_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466.vcf.gz").getAbsolutePath();
            String outfileS = new File(outfileDirS,"chr" + chrArr[j] + "_vmap2.1_heter_SNPbased_DomesticatedEmmer_PI355466_chrposGenotype.txt.gz").getAbsolutePath();
            String logfileS = new File(logDirS,"log_" + new File(outfileS).getName().split(".gz")[0]).getAbsolutePath();

            System.out.println("nohup java -jar 034_mkindividualVCFtoChrPosGenotype.jar " + infileS + " " + outfileS + " > " + logfileS  + " 2>&1 &" );
        }

        //AABBDD 一个cultivar个体的运行
//        nohup java -jar 034_mkindividualVCFtoChrPosGenotype.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/003_indiVCF/chr1B_vmap2.1_heter_SNPbased_Cultivar.vcf.gz /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr1B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/003/log_chr1B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt 2>&1 &
//        nohup java -jar 034_mkindividualVCFtoChrPosGenotype.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/003_indiVCF/chr2A_vmap2.1_heter_SNPbased_Cultivar.vcf.gz /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr2A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/003/log_chr2A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt 2>&1 &
//        nohup java -jar 034_mkindividualVCFtoChrPosGenotype.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/003_indiVCF/chr2B_vmap2.1_heter_SNPbased_Cultivar.vcf.gz /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr2B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/003/log_chr2B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt 2>&1 &
//        nohup java -jar 034_mkindividualVCFtoChrPosGenotype.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/003_indiVCF/chr3B_vmap2.1_heter_SNPbased_Cultivar.vcf.gz /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr3B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/003/log_chr3B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt 2>&1 &
//        nohup java -jar 034_mkindividualVCFtoChrPosGenotype.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/003_indiVCF/chr6B_vmap2.1_heter_SNPbased_Cultivar.vcf.gz /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr6B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/003/log_chr6B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt 2>&1 &
//        nohup java -jar 034_mkindividualVCFtoChrPosGenotype.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/003_indiVCF/chr7A_vmap2.1_heter_SNPbased_Cultivar.vcf.gz /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr7A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/003/log_chr7A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt 2>&1 &
//        nohup java -jar 034_mkindividualVCFtoChrPosGenotype.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/003_indiVCF/chr7B_vmap2.1_heter_SNPbased_Cultivar.vcf.gz /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr7B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/003/log_chr7B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt 2>&1 &
        //
        //
        //nohup java -jar 034_mkindividualVCFtoChrPosGenotype.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/003_indiVCF/chr1A_vmap2.1_heter_SNPbased_Cultivar.vcf.gz /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr1A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/003/log_chr1A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt 2>&1 &
        //nohup java -jar 034_mkindividualVCFtoChrPosGenotype.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/003_indiVCF/chr1D_vmap2.1_heter_SNPbased_Cultivar.vcf.gz /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr1D_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/003/log_chr1D_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt 2>&1 &
        //nohup java -jar 034_mkindividualVCFtoChrPosGenotype.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/003_indiVCF/chr2D_vmap2.1_heter_SNPbased_Cultivar.vcf.gz /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr2D_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/003/log_chr2D_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt 2>&1 &
        //nohup java -jar 034_mkindividualVCFtoChrPosGenotype.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/003_indiVCF/chr3A_vmap2.1_heter_SNPbased_Cultivar.vcf.gz /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr3A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/003/log_chr3A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt 2>&1 &
        //nohup java -jar 034_mkindividualVCFtoChrPosGenotype.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/003_indiVCF/chr3D_vmap2.1_heter_SNPbased_Cultivar.vcf.gz /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr3D_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/003/log_chr3D_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt 2>&1 &
        //nohup java -jar 034_mkindividualVCFtoChrPosGenotype.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/003_indiVCF/chr4A_vmap2.1_heter_SNPbased_Cultivar.vcf.gz /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr4A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/003/log_chr4A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt 2>&1 &
        //nohup java -jar 034_mkindividualVCFtoChrPosGenotype.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/003_indiVCF/chr4B_vmap2.1_heter_SNPbased_Cultivar.vcf.gz /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr4B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/003/log_chr4B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt 2>&1 &
        //nohup java -jar 034_mkindividualVCFtoChrPosGenotype.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/003_indiVCF/chr4D_vmap2.1_heter_SNPbased_Cultivar.vcf.gz /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr4D_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/003/log_chr4D_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt 2>&1 &
        //nohup java -jar 034_mkindividualVCFtoChrPosGenotype.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/003_indiVCF/chr5A_vmap2.1_heter_SNPbased_Cultivar.vcf.gz /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr5A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/003/log_chr5A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt 2>&1 &
        //nohup java -jar 034_mkindividualVCFtoChrPosGenotype.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/003_indiVCF/chr5B_vmap2.1_heter_SNPbased_Cultivar.vcf.gz /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr5B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/003/log_chr5B_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt 2>&1 &
        //nohup java -jar 034_mkindividualVCFtoChrPosGenotype.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/003_indiVCF/chr5D_vmap2.1_heter_SNPbased_Cultivar.vcf.gz /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr5D_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/003/log_chr5D_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt 2>&1 &
        //nohup java -jar 034_mkindividualVCFtoChrPosGenotype.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/003_indiVCF/chr6A_vmap2.1_heter_SNPbased_Cultivar.vcf.gz /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr6A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/003/log_chr6A_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt 2>&1 &
        //nohup java -jar 034_mkindividualVCFtoChrPosGenotype.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/003_indiVCF/chr6D_vmap2.1_heter_SNPbased_Cultivar.vcf.gz /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr6D_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/003/log_chr6D_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt 2>&1 &
        //nohup java -jar 034_mkindividualVCFtoChrPosGenotype.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/003_indiVCF/chr7D_vmap2.1_heter_SNPbased_Cultivar.vcf.gz /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr7D_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/003/log_chr7D_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt 2>&1 &

// DD 一个个体的运行
        //nohup java -jar 034_mkindividualVCFtoChrPosGenotype.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/003_indiVCF/chr1D_vmap2.1_heter_SNPbased_Ae.tauschii_AE1211.vcf.gz /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr1D_vmap2.1_heter_SNPbased_Ae.tauschii_AE1211_chrposGenotype.txt.gz > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/003/log_chr1D_vmap2.1_heter_SNPbased_Ae.tauschii_AE1211_chrposGenotype.txt 2>&1 &
        //nohup java -jar 034_mkindividualVCFtoChrPosGenotype.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/003_indiVCF/chr2D_vmap2.1_heter_SNPbased_Ae.tauschii_AE1211.vcf.gz /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr2D_vmap2.1_heter_SNPbased_Ae.tauschii_AE1211_chrposGenotype.txt.gz > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/003/log_chr2D_vmap2.1_heter_SNPbased_Ae.tauschii_AE1211_chrposGenotype.txt 2>&1 &
        //nohup java -jar 034_mkindividualVCFtoChrPosGenotype.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/003_indiVCF/chr3D_vmap2.1_heter_SNPbased_Ae.tauschii_AE1211.vcf.gz /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr3D_vmap2.1_heter_SNPbased_Ae.tauschii_AE1211_chrposGenotype.txt.gz > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/003/log_chr3D_vmap2.1_heter_SNPbased_Ae.tauschii_AE1211_chrposGenotype.txt 2>&1 &
        //nohup java -jar 034_mkindividualVCFtoChrPosGenotype.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/003_indiVCF/chr4D_vmap2.1_heter_SNPbased_Ae.tauschii_AE1211.vcf.gz /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr4D_vmap2.1_heter_SNPbased_Ae.tauschii_AE1211_chrposGenotype.txt.gz > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/003/log_chr4D_vmap2.1_heter_SNPbased_Ae.tauschii_AE1211_chrposGenotype.txt 2>&1 &
        //nohup java -jar 034_mkindividualVCFtoChrPosGenotype.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/003_indiVCF/chr5D_vmap2.1_heter_SNPbased_Ae.tauschii_AE1211.vcf.gz /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr5D_vmap2.1_heter_SNPbased_Ae.tauschii_AE1211_chrposGenotype.txt.gz > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/003/log_chr5D_vmap2.1_heter_SNPbased_Ae.tauschii_AE1211_chrposGenotype.txt 2>&1 &
        //nohup java -jar 034_mkindividualVCFtoChrPosGenotype.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/003_indiVCF/chr6D_vmap2.1_heter_SNPbased_Ae.tauschii_AE1211.vcf.gz /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr6D_vmap2.1_heter_SNPbased_Ae.tauschii_AE1211_chrposGenotype.txt.gz > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/003/log_chr6D_vmap2.1_heter_SNPbased_Ae.tauschii_AE1211_chrposGenotype.txt 2>&1 &
        //nohup java -jar 034_mkindividualVCFtoChrPosGenotype.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/003_indiVCF/chr7D_vmap2.1_heter_SNPbased_Ae.tauschii_AE1211.vcf.gz /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/005_003_trans/chr7D_vmap2.1_heter_SNPbased_Ae.tauschii_AE1211_chrposGenotype.txt.gz > /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/003/log_chr7D_vmap2.1_heter_SNPbased_Ae.tauschii_AE1211_chrposGenotype.txt 2>&1 &


    }



    /**
     * 将单个样品的VCF文件转化为可以计算片段杂合度的文件格式；
     * 0/0 为0； 0/1为1； 1/1 为2； ./. 为NA
     *
     */
    public void mkGenotype(String infileS, String outfileS){
//        String infileS= "/Users/Aoyue/Documents/0/chr1A_vmap2.1_heter_SNPbased_Cultivar.vcf.gz";
//        String outfileS = "/Users/Aoyue/Documents/1/chr1A_vmap2.1_heter_SNPbased_Cultivar_heter.txt.gz";

        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            bw.write("Chr\tPos\tGenotype");
            bw.newLine();


            //特殊情况，VCFTOOL把log文件也读到文件中去了，所以要过滤开头那几行
//            for (int i = 0; i < 13 ; i++) {
//                br.readLine();
//            }
            int cnt =0;
            List<String> l = new ArrayList<>();
            String temp = null;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) {
                } else {
                    l= PStringUtils.fastSplit(temp);
                    String chr = l.get(0);
                    String pos = l.get(1);
                    String genoArray = l.get(9);
                    String geno = PStringUtils.fastSplit(genoArray,":").get(0);
//                    System.out.println(geno);
                    if(geno.equals("0/0")){
                        bw.write(chr + "\t" + pos + "\t0");
                        bw.newLine();
                        cnt++;
                    }
                    if(geno.equals("0/1")){
                        bw.write(chr + "\t" + pos + "\t1");
                        bw.newLine();
                        cnt++;
                    }
                    if(geno.equals("1/1")){
                        bw.write(chr + "\t" + pos + "\t2");
                        bw.newLine();
                        cnt++;
                    }
                    if(geno.equals("./.")){
                        bw.write(chr + "\t" + pos + "\tNA");
                        bw.newLine();
                        cnt++;
                    }

                }
            }
            bw.flush();
            bw.close();
            br.close();
            System.out.println(infileS + " is completed at " + outfileS);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }


    public void mergeTxt(){
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/004_heterogozysity/003_indi_test/003_cal2MWindow_1Mstep_2";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/004_heterogozysity/003_indi_test/004_merge2/heter_cultivar_indi_2Mwindow_1Mstep.txt";
//        new CountSites().mergeTxt(infileDirS,outfileS);

    }


    /**
     * 根据点的数值，计算 滑窗的数值
     */
    public void script_calWindowStep(){

        //方法1：循环法
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/004_heterogozysity/003_indi_test/002_out2";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/004_heterogozysity/003_indi_test/003_cal2MWindow_1Mstep_2";
//        String[] chrArr = {"1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B", "7D"};
//        for (int i = 0; i < chrArr.length; i++) {
////            String infileS = new File(infileDirS,"chr" + chrArr[i] + "_vmap2.1_heter_SNPbased_Cultivar.txt.gz").getAbsolutePath();
////            String outfileS = new File(outfileDirS,"chr" + chrArr[i] + "_vmap2.1_heter_SNPbased_Cultivar_2Mwindow_1Mstep.txt").getAbsolutePath();
//
//            String infileS = new File(infileDirS,"chr" + chrArr[i] + "_vmap2.1_heter_SNPbased_Landrace.txt.gz").getAbsolutePath();
//            String outfileS = new File(outfileDirS,"chr" + chrArr[i] + "_vmap2.1_heter_SNPbased_Landrace_2Mwindow_1Mstep.txt").getAbsolutePath();
//
//            HashMap<Integer,String> hm = new AoFile().getHashMapintKey(infileS,1,2);
//            new Bin().calwindowstep(chrArr[i],hm,2000000,1000000,outfileS);
//        }


        //方法2：列表法

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/004_heterogozysity/003_indi_test/002_out2";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/004_heterogozysity/003_indi_test/003_cal2MWindow_1Mstep_2";
        List<File> fsList = IOUtils.getVisibleFileListInDir(infileDirS);
        fsList.parallelStream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                String outfileS = null;
                BufferedReader br = null;
                if (infileS.endsWith(".txt")) {
                    br = IOUtils.getTextReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().split(".txt")[0] + "_2Mwindow_1Mstep.txt").getAbsolutePath();
                } else if (infileS.endsWith(".txt.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().split(".txt.gz")[0] + "_2Mwindow_1Mstep.txt").getAbsolutePath();
                }

                String chr = f.getName().substring(3,5);
                HashMap<Integer,String> hm = new AoFile().getHashMapintKey(infileS,1,2);
                new Bin().calwindowstep(chr,hm,2000000,1000000,outfileS);

                System.out.println(f.getName() + "\tis completed at " + outfileS);
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }


    /**
     * 提取个体的0/1位点信息，以及杂合度标记为1
     * 这里只是把一个个体的基因型是 0/1 的位点找出来了，并没有做其他任何措施
     */
    public void scriptforIndi(){
        //程序运行时，输入输出路径设置
        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/013_VMapIIbyRef";
        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/002_out";
//        String taxafileS = "/data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/001_taxalist/Cultivar.txt";
        String taxafileS = "/data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/004_out_indivi/001_taxalist/Landrace.txt";

//        String group = "Cultivar";
        String group = "Landrace";

        String logDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/002";

        String[] chrArr = {"1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B", "7D"};
        for (int j = 0; j < chrArr.length; j++) {
            String infileS = new File(infileDirS,"chr" + chrArr[j] + "_vmap2.1.vcf").getAbsolutePath();
            String outfileS = new File(outfileDirS,"chr" + chrArr[j] + "_vmap2.1_heter_SNPbased_" + group + ".txt.gz").getAbsolutePath();
            String logfileS = new File(logDirS,"log_chr" + chrArr[j] + "_vmap2.1_heter_SNPbased_" + group + ".txt").getAbsolutePath();
            System.out.println("java -jar 033_getSNPHeterbySite.jar " + infileS + " " + outfileS + " " + taxafileS + " > " + logfileS  + " &" );
        }

//        String[] chrArr = {"1A", "1B", "2A", "2B", "3A", "3B", "4A", "4B", "5A", "5B", "6A", "6B", "7A", "7B"};
//        for (int j = 0; j < chrArr.length; j++) {
//            String infileS = new File(infileDirS,"chr" + chrArr[j] + "_vmap2.1.vcf").getAbsolutePath();
//            String outfileS = new File(outfileDirS,"chr" + chrArr[j] + "_vmap2.1_heter_SNPbased_" + group + ".txt.gz").getAbsolutePath();
//            String logfileS = new File(logDirS,"log_chr" + chrArr[j] + "_vmap2.1_heter_SNPbased_" + group + ".txt").getAbsolutePath();
//            System.out.println("java -jar 033_getSNPHeterbySite.jar " + infileS + " " + outfileS + " " + taxafileS + " > " + logfileS  + " &" );
//        }

//        String[] chrArr = {"1D","2D", "3D", "4D", "5D", "6D","7D"};
//        for (int j = 0; j < chrArr.length; j++) {
//            String infileS = new File(infileDirS,"chr" + chrArr[j] + "_vmap2.1.vcf").getAbsolutePath();
//            String outfileS = new File(outfileDirS,"chr" + chrArr[j] + "_vmap2.1_heter_SNPbased_" + group + ".txt.gz").getAbsolutePath();
//            String logfileS = new File(logDirS,"log_chr" + chrArr[j] + "_vmap2.1_heter_SNPbased_" + group + ".txt").getAbsolutePath();
//            System.out.println("java -jar 033_getSNPHeterbySite.jar " + infileS + " " + outfileS + " " + taxafileS + " > " + logfileS  + " &" );
//        }


    }

    /**
     * 进行window step的小测试
     */
    public void windowCal(){
//        String infileS= "/Users/Aoyue/Documents/chr002.subgenome.maf0.01.SNP_bi.cultivar.vcf.txt";
//        String outfileS = "/Users/Aoyue/Documents/chr002_Cultivar_100kwindow50kstep.txt";
//        HashMap<Integer,String> hm = new AoFile().getHashMapintKey(infileS,1,2);
        //        new Bin().calwindow("2",hm,1000000,outfileS);
//        new Bin().calwindowstep("2",hm,100000,50000,outfileS);


        String infileS= "/Users/Aoyue/Documents/a.txt";
        String outfileS = "/Users/Aoyue/Documents/b.txt";
        HashMap<Integer,String> hm = new AoFile().getHashMapintKey(infileS,1,2);
//        new Bin().calwindowstep_ResidualHeterozygosity("2A",hm,2000000,1000000,outfileS);

    }

    /**
     * 计算不同倍性的小麦的杂合子，只保留没有分离的位点。
     */
    public void scriptSNPbased(){

        //***************************** step one : 确定其倍性，根据倍性计算 ****************************//
        String dbfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/000_prepareData/001_input/taxaList.txt";
        new AoFile().readheader(dbfileS);
        HashMap<String,String> hm = new AoFile().getHashMapStringKey(dbfileS,10,8);
        System.out.println(hm.entrySet());
        List<String> groupl = new ArrayList<String>(hm.keySet());
        Collections.sort(groupl);
        //***************************** step two  ****************************//

        //程序运行时，输入输出路径设置
        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/013_VMapIIbyRef";
        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/003_out";
        String taxaDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/001_treeValidatedGroup_byPloidy";
        String logDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/log/001";

        String scriptS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/004_heterogozysity/002_script_SNPbased/sh_heter20200123.sh";

        //java -jar 033_getSNPHeterbySite.jar /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/test/chr002.subgenome.maf0.01.SNP_bi.subset.vcf /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/test/chr002.subgenome.maf0.01.SNP_bi.cultivar.heter.txt /data4/home/aoyue/vmap2/analysis/021_popGen/004_heter/test/Cultivar.txt &;
        //chr1A_vmap2.1.vcf

        try {

            BufferedWriter bw = IOUtils.getTextWriter(scriptS);
            for (int i = 0; i < groupl.size(); i++) {
                String group = groupl.get(i);
                if(group.equals("ExclusionHexaploid") || group.equals("ExclusionTetraploid"))continue;
                if(group.equals("Hexaploid")){
                    String taxafileS = new File(taxaDirS,group+".txt").getAbsolutePath();
                    String[] chrArr = {"1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B", "7D"};
                    for (int j = 0; j < chrArr.length; j++) {
                        String infileS = new File(infileDirS,"chr" + chrArr[j] + "_vmap2.1.vcf").getAbsolutePath();
                        String outfileS = new File(outfileDirS,"chr" + chrArr[j] + "_vmap2.1_heter_SNPbased_" + group + ".txt.gz").getAbsolutePath();
                        String logfileS = new File(logDirS,"log_chr" + chrArr[j] + "_vmap2.1_heter_SNPbased_" + group + ".txt").getAbsolutePath();
                        System.out.println("java -jar 033_getSNPHeterbySite.jar " + infileS + " " + outfileS + " " + taxafileS + " > " + logfileS );
                        bw.write("java -jar 033_getSNPHeterbySite.jar " + infileS + " " + outfileS + " " + taxafileS + " > " + logfileS);
                        bw.newLine();
                    }
                }
                if(group.equals("Tetraploid")){
                    String[] chrArr = {"1A", "1B", "2A", "2B", "3A", "3B", "4A", "4B", "5A", "5B", "6A", "6B", "7A", "7B"};
                    String taxafileS = new File(taxaDirS,group+".txt").getAbsolutePath();
                    for (int j = 0; j < chrArr.length; j++) {
                        String infileS = new File(infileDirS,"chr" + chrArr[j] + "_vmap2.1.vcf").getAbsolutePath();
                        String outfileS = new File(outfileDirS,"chr" + chrArr[j] + "_vmap2.1_heter_SNPbased_" + group + ".txt.gz").getAbsolutePath();
                        String logfileS = new File(logDirS,"log_chr" + chrArr[j] +  "_vmap2.1_heter_SNPbased_" + group + ".txt").getAbsolutePath();
                        System.out.println("java -jar 033_getSNPHeterbySite.jar " + infileS + " " + outfileS + " " + taxafileS + " > " + logfileS );
                        bw.write("java -jar 033_getSNPHeterbySite.jar " + infileS + " " + outfileS + " " + taxafileS + " > " + logfileS);
                        bw.newLine();
                    }


                }
                if(group.equals("Ae.tauschii")){
                    String[] chrArr = {"1D","2D", "3D", "4D", "5D", "6D","7D"};
                    String taxafileS = new File(taxaDirS,group+".txt").getAbsolutePath();
                    for (int j = 0; j < chrArr.length; j++) {
                        String infileS = new File(infileDirS,"chr" + chrArr[j] + "_vmap2.1.vcf").getAbsolutePath();
                        String outfileS = new File(outfileDirS,"chr" + chrArr[j] + "_vmap2.1_heter_SNPbased_" + group + ".txt.gz").getAbsolutePath();
                        String logfileS = new File(logDirS,"log_chr" + chrArr[j] + "_vmap2.1_heter_SNPbased_" + group + ".txt").getAbsolutePath();
                        System.out.println("java -jar 033_getSNPHeterbySite.jar " + infileS + " " + outfileS + " " + taxafileS + " > " + logfileS );
                        bw.write("java -jar 033_getSNPHeterbySite.jar " + infileS + " " + outfileS + " " + taxafileS + " > " + logfileS);
                        bw.newLine();
                    }
                }
            }

            bw.flush();
            bw.close();

            new SplitScript().splitScript(scriptS,"heter",14,3);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
    
}


