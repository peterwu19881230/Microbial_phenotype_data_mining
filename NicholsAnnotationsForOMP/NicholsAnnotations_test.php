<?php
/*
This code will be used to insert annotations into the website after strain pages are created.

Execute by: php NicholsAnnotations_test.php -w /var/www/html/omp/sandy

The test website's URL: https://microbialphenotypes.org//sandy/index.php?title=Main_Page


*/
$params = getopt( "w:" );
set_IP($params['w']);
$maintClass = "DemoMaint";


class DemoMaint extends Maintenance {
    public function __construct() {
        parent::__construct();
        $this->parse_parameters();
    }
    #I start by opening the file that contains all of the S. pombe annotations
    public function execute() {


$parentStrain="OMP_ST:800_!_Escherichia_coli_K-12_BW25113";


            #Get the strain page
            $dbr = wfGetDB( DB_SLAVE );
            $result = $dbr->select(
                array('page'),
                '*',
                array("page_title LIKE '$parentStrain'"), #In MySQL " " is "_"  #Have to have "-" otherwise GAL10 and GAL1 are found at the same time
                __METHOD__,
                array()
            );



            foreach($result as $x){

                #Uses the pages found and creates a title object and finds the tables on the strain pages that were previously created using the pombe_knockout_insert.php code
                $newtitle = Title::newFromText($x->page_title);
                $new_page = new WikiPageTE($newtitle);
                echo $newtitle."\n";


                #Make an annotation based on the page found
                $strain_table_template = "OMP_annotation_table";
                $newtable = $new_page->getTable($strain_table_template);
                $box = $newtable[0];
                $newrow = $box->insert_row(
                          # Annotation ID: auto-generated
                    "||". # Qualifier
                    "||". "OMP:0007731". #OMP ID
                    "||"."abolished vegetative cell population viability". #OMP term name //For this guy we have to change the term (not creating) to "population" something
                    "||"."OMP_AN:24800". #Relative phenotype information: get annotation id from parent strain and put in relative phenotype box #Relative phenotype information: Relative to. #The following Genotype differences and Condition differences are auto-generated
                    "||". "medium:"."YPD".#Experimental condition
                    "||"."ECO:0005004".#ECO ID
                    "||"."cell viability assay evidence ". #ECO term name
                    "||"."PMID:"."123456". #Reference: Pubmed ID
                    "||".#Annotation Extension: put the CHEBI id here
                    "||" #Notes
                );
                $new_page->touch();
                echo "An annotation has been made\n";



            }


    }
    private function parse_parameters(){
        $this->addOption( "people", "say hello to this", $required = false, $withArg = true, $shortName = 'l' );
    }
}
require_once( RUN_MAINTENANCE_IF_MAIN );
function set_IP($path){
    global $IP;
    if ( isset($path) && is_file("$path/maintenance/Maintenance.php") ){
        $IP = $path;
        require_once( $IP . "/maintenance/Maintenance.php" );
        return $path;
    } else {
        die ("need -w <path-to-wiki-directory>");
    }
}				