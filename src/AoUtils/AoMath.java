/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package AoUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import utils.IOUtils;
import utils.PStringUtils;

/**
 *
 * @author Aoyue
 */
public class AoMath {

    public AoMath() {
        this.filterValue();

    }

    /**
     * 过滤DAF_ABD等于0或者1的位点
     *
     */
    public void filterValue() {

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/004_merge";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/005_filterDAF_ABD";

        //  new CountSites().mergesubsetVCF(args[0], args[1]);
        // new CountSites().mergesubsetVCF(args[0], args[1], args[2]);
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                System.out.println(fs[i].getName() + " is hidden");
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        fsList.stream().forEach(f -> {

            try {
                String infileS = f.getAbsolutePath();
                String outfileS = new File(outfileDirS, f.getName().replaceFirst(".txt", "_filterDAF_ABD10.txt")).getAbsolutePath();
                BufferedReader br = IOUtils.getTextReader(infileS);
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                String temp = null;
                String header = br.readLine();
                bw.write(header);
                bw.newLine();
                int cnt = 0;
                List<String> l = new ArrayList();
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    StringBuilder sb = new StringBuilder();
                    cnt++;
                    String DAF_ABD = l.get(15);
                    if (DAF_ABD.startsWith("N")) {
                        continue;
                    }
                    double daf_abd = Double.parseDouble(DAF_ABD);
                    if (daf_abd == 0 || (daf_abd == 1)) {
                        continue;
                    }
                    bw.write(temp);
                    bw.newLine();

                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println();
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }

        });
        /*==================================== 测试用 =============================================*/

    }

}
