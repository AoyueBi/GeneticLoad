/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Plot;

import format.table.RowTable;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.HashMap;
import java.util.List;
import utils.IOUtils;
import utils.PStringUtils;

/**
 *
 * @author Aoyue
 */
public class PCA {
    
    public PCA(){
        
    }
    
    
    public void addGrouptoMDS_AandBsubgenome(){
        
        
    }
            
    
    // 对MDS方法生成的PC结果添加分组信息，在第二列中加入 
    //String[] group = {"Oceania","Africa","North America","South America","Europe","Central Asia","South Asia","Western Asia","East Asia","NA"};
    //String[] col = {"#F1E1FF","#F4D03F","#F1948A","#5DADE2","#ABEBC6","#239B56","#CD6155","#FF6347","#7B241C","#EBEDEF"};
    //String[] col = {"#5DADE2","#7B241C","#F1E1FF","#F4D03F","#FF9900","#006600","#389038","#82C782","#CCFFCC","#EBEDEF"};
    //分组1为：indexColum 
    //分组2为：大洋洲 非洲 北美洲 南美洲 欧洲 亚洲 
    //"Oceania","Africa","North America","South America","Europe","Asia" 
    //颜色为："#F1E1FF","#F4D03F","#5DADE2","#7B241C","#FF9900","#82C782"
    // 分组3为：大洋洲 非洲 北美洲 南美洲 欧洲部洲 亚洲部洲
    //分组4为国家
    //0 <- "Africa" 
    // 1 <- "Asia" 
    // 2 <- "Europe" 
    // 3 <- "North America" 
    // 4 <- * "Oceania"
    // 5 <- "South America"
    //6<- "NA"
    public void addGrouptoMDS() {
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/005_pca/ab/001_mdsMethod/002_MDS_PCs_Matrix_subset63ksnp_forR.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/005_pca/ab/001_mdsMethod/003_MDS_PCs_Matrix_subset63ksnp_forR_addGroup.txt";
        String dbfileS = "/Users/Aoyue/project/wheatVMapII/001_germplasm/GermplasmDB/wheatVMapII_AB_S205_germplasmInfo.txt";
        HashMap<String, String> hmDatabaseIDContinent = new HashMap<>();
        HashMap<String, String> hmDatabaseIDPartContinent = new HashMap<>();
        HashMap<String, String> hmDatabaseIDCty = new HashMap<>();
        HashMap<String, Integer> hmContinentIndex = new HashMap<>();
        HashMap<String, String> hmEmmertype = new HashMap<>();

        RowTable<String> t = new RowTable<>(dbfileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            String databaseID = t.getCell(i, 0);
            String continent = t.getCell(i, 12);
            String partContinent = t.getCell(i, 14);
            String country = t.getCell(i, 10);
            String type = t.getCell(i, 7);
            hmDatabaseIDContinent.put(databaseID, continent);
            hmDatabaseIDPartContinent.put(databaseID, partContinent);
            hmDatabaseIDCty.put(databaseID, country);
            hmEmmertype.put(databaseID, type);
        }

        String[] continents = {"Africa", "Asia", "Europe", "North America", "Oceania", "South America"};
        for (int i = 0; i < continents.length; i++) {
            hmContinentIndex.put(continents[i], i);
        }

        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = br.readLine(); //header
            List<String> l = PStringUtils.fastSplit(temp);
            bw.write(l.get(0) + "\t");
            bw.write("GroupIndex\tCountry\tContinent\tPart-Cpntinent\tType");
            for (int i = 1; i < l.size(); i++) {
                bw.write("\t" + l.get(i));
            }
            bw.newLine();
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String id = l.get(0);
                String cty = hmDatabaseIDCty.get(id);
                if (id.equals("CS")) {
                    bw.write(l.get(0) + "\t");
                    bw.write("1\tCHN\tAsia\tEast Asia\tdicoccum");
                    for (int i = 1; i < l.size(); i++) {
                        bw.write("\t" + l.get(i));
                    }
                    bw.newLine();

                } else {
                    if (cty.equals("NA")) {
                        bw.write(l.get(0) + "\t");
                        bw.write("6\tNA\tNA\tNA\t" + hmEmmertype.get(l.get(0)));
                        for (int i = 1; i < l.size(); i++) {
                            bw.write("\t" + l.get(i));
                        }
                        bw.newLine();

                    } else {
                        bw.write(id + "\t");
                        bw.write(String.valueOf(hmContinentIndex.get(hmDatabaseIDContinent.get(id))) + "\t" + hmDatabaseIDCty.get(id) + "\t" + hmDatabaseIDContinent.get(id)
                                + "\t" + hmDatabaseIDPartContinent.get(id) + "\t" + hmEmmertype.get(l.get(0)));
                        for (int i = 1; i < l.size(); i++) {
                            bw.write("\t" + l.get(i));
                        }
                        bw.newLine();

                    }

                }

            }
            br.close();
            bw.flush();
            bw.close();

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
    
}
