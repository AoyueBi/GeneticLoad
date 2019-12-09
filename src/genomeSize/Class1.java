/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.wheat.genomeSize;

/**
 *
 * @author mingzhang
 */
public class Class1 {
    protected InnerClass1 ic;
    
    public Class1(){
        ic = new InnerClass1();
}
    public void displayString(){
        System.out.println(ic.getString() + ".");
        System.out.println(ic.getAnotherString() + ".");
    }
    
    static public void main(String[] args){
       Class1 c1 = new Class1();
       c1.displayString();
    }
    
    protected class InnerClass1{
        public String getString(){
            return "InnerClass1:getString invoked";
        }
        
        public String getAnotherString(){
            return "InnerClass1: getAnotherString invoked";
        }
    }
    
    
}
