package PopulationAnalysis;

import AoUtils.SplitScript;
import pgl.utils.IOUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

public class Fst {
    public Fst() {
//        this.mkFstCommandbasedSNP();
        new SplitScript().splitScript2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/101_Fst/001_scriptSNPbased/fst_basedSNP_20200205.sh",7,23);

    }

    /**
     *
     */
    public void mkFstCommandbasedSNP() {

        //local file： one: group two: script
        String groupHexaandTetraDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/101_Fst/000_group/hexaandTetra";
        String groupHexaandDiDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/101_Fst/000_group/hexaandDi";
        String shScriptDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/101_Fst/000_scriptSNPbased";

        // HPC file: group fileDirS
        String group1FileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/101_Fst/000_group/hexaandTetra";
        String group2FileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/101_Fst/000_group/hexaandDi";

        // HPC file: output fileDirS
        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/013_VMapIIbyRef";
        String outputDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/101_Fst/002_fst_basedSNP";

        List<File> fs = IOUtils.getVisibleFileListInDir(groupHexaandTetraDirS);
        File[] group1FileS = fs.toArray(new File[fs.size()]);

        List<File> fs2 = IOUtils.getVisibleFileListInDir(groupHexaandDiDirS);
        File[] group2FileS = fs2.toArray(new File[fs2.size()]);

        new File(shScriptDirS).mkdirs(); //无论在不在，都可以创建该文件夹

        ArrayList<String> perlList = new ArrayList(); //在循环外建立perlList集合， 每个集合包含多个字符串，一个字符串代表一个文件。

        for (int i = 0; i < group1FileS.length - 1; i++) {
            String pop1 = group1FileS[i].getName().replace(".txt", ""); //第一组的名字
            for (int j = i + 1; j < group1FileS.length; j++) {
                String pop2 = group1FileS[j].getName().replace(".txt", ""); //第二组的名字
                //vcftools --gzvcf test.vcf.gz --weir-fst-pop ../groups/Teosinte.txt --weir-fst-pop ../groups/Stiff_stalk.txt --weir-fst-pop ../groups/Non_stiff_stalk.txt --out out.txt
                String[] chrArr = {"1A", "1B", "2A", "2B", "3A", "3B", "4A", "4B", "5A", "5B", "6A", "6B", "7A", "7B"};
                for (int k = 0; k < chrArr.length; k++) {
                    String chr = chrArr[k];
                    String infileS = new File(infileDirS, "chr" + chr + "_vmap2.1.vcf").getAbsolutePath();
                    String outfileS = new File(outputDirS, pop1 + "_VS_" + pop2 + "_chr" + chr).getAbsolutePath();
                    String group1S = new File(group1FileDirS, group1FileS[i].getName()).getAbsolutePath();
                    String group2S = new File(group1FileDirS, group1FileS[j].getName()).getAbsolutePath();
                    System.out.println("vcftools --vcf " + infileS + " --weir-fst-pop " + group1S + " --weir-fst-pop " + group2S + " --out " + outfileS);
                }
            }
        }

        for (int i = 0; i < group2FileS.length - 1; i++) {
            String pop1 = group2FileS[i].getName().replace(".txt", "");
            for (int j = i + 1; j < group2FileS.length; j++) {
                String pop2 = group2FileS[j].getName().replace(".txt", "");
                String[] chrArr = {"1D", "2D", "3D", "4D", "5D", "6D", "7D"};
                for (int k = 0; k < chrArr.length; k++) {
                    String chr = chrArr[k];
                    String infileS = new File(infileDirS, "chr" + chr + "_vmap2.1.vcf").getAbsolutePath();
                    String outfileS = new File(outputDirS, pop1 + "_VS_" + pop2 + "_chr" + chr).getAbsolutePath();
                    String group1S = new File(group2FileDirS, group2FileS[i].getName()).getAbsolutePath();
                    String group2S = new File(group2FileDirS, group2FileS[j].getName()).getAbsolutePath();
                    System.out.println("vcftools --vcf " + infileS + " --weir-fst-pop " + group1S + " --weir-fst-pop " + group2S + " --out " + outfileS);
                }
            }
        }
    }

}



