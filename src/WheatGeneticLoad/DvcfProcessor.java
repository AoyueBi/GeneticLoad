/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package WheatGeneticLoad;

import AoUtils.CountSites;

/**
 *
 * @author Aoyue
 */
public class DvcfProcessor {


    public DvcfProcessor() {
        this.step1();
        
    }
    
    public void step1(){
        CountSites c = new CountSites();
//        c.countSites_parallelStream("/data1/home/aoyue/fastcall_Dgenome/rawVCF");  //java -jar PlantGenetics.jar > countSites_Dgenome.txt &  //deprecated
//        c.mergeChr1and2_Dgenome("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/D/001_countSites/countSites_Dgenome.txt", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/D/001_countSites/countSites_mergeChr1and2_Dgenome.txt"); //deprecated
//        c.countSitesinFastCallformat("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/D/test");
//        c.filterSNPtoBi("/data4/home/aoyue/vmap2/analysis/001_rawvcf/d/", "/data4/home/aoyue/vmap2/analysis/002_bivcf/d/");
        
        //c.filterSNPtoBi("/data1/home/aoyue/fastcall_Dgenome/rawVCF/", "/data1/home/aoyue/fastcall_Dgenome/002_bivcf/d/");
        //c.subsetVCFRandomParallel_GZ("/data4/home/aoyue/vmap2/analysis/003_filterMiss/d/", "/data4/home/aoyue/vmap2/analysis/004_subsetvcf/d/");
        //c.subsetVCFRandomParallel("/Users/Aoyue/Documents/test/", "/Users/Aoyue/Documents/out/");
    }
}
