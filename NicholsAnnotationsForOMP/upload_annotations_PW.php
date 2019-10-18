<?php
/*
This code will be used to insert annotations into the website after strain pages are created.

Execute by: php NicholsStrainsCreation.php -w /Library/WebServer/Documents/omp/peter/ -f file_for_annotations
*/


$params = getopt( "w:" ); #-w to point to the path of the wiki. For other parameters, don't do it here. Do it in parse_parameters()
set_IP($params['w']);

$maintClass = "DemoMaint";
class DemoMaint extends Maintenance {
    public function __construct() {
        parent::__construct();
        $this->parse_parameters();
    }


    public function execute() {
        #this opens the pombase file
        $filename = $this->getOption('file_name');
        echo "opening $filename\n";
        $fh = fopen($filename, 'r');
        if(!$fh) die;



        for ($i=0;;$i++){


            $line = fgets($fh);

            if ($line{0} == "\n" ||trim($line) == ""){
                echo("---- Annotation completed ----\n");
                break;
            }




            #this makes variables for all of the data from the file
            list ($page_title,$OMP_ID,$OMP_term_name,$relative_phenotype_info,$Exp_condition,$ECO_ID,$ECO_term_name,$PMID,$Annot_extension,$Notes)= explode("\t",trim($line));



            #once I have the paper I am looking for, I search by allele, which in this case is also the strain name
            $dbr = wfGetDB( DB_SLAVE );
            $result = $dbr->select(
                array('page'),
                '*',
                array("page_title = '$page_title' "),
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
                    "||". $OMP_ID. #OMP ID
                    "||".$OMP_term_name. #OMP term name //For this guy we have to change the term (not creating) to "population" something
                    "||".$relative_phenotype_info. #Relative phenotype information: get annotation id from parent strain and put in relative phenotype box #Relative phenotype information: Relative to. #The following Genotype differences and Condition differences are auto-generated
                    "||".$Exp_condition.#Experimental condition
                    "||".$ECO_ID.#ECO ID
                    "||".$ECO_term_name. #ECO term name
                    "||".$PMID. #Reference: Pubmed ID
                    "||".$Annot_extension. #Annotation Extension: put the CHEBI id here
                    "||".$Notes #Notes
                );
                $new_page->touch();
                echo "An annotation has been made\n";



            }



        }
    }
    private function parse_parameters(){ #Ref: https://www.mediawiki.org/wiki/Manual:Writing_maintenance_scripts
        $this->addOption( "file_name", "put your file name", $required = true, $withArg = true, $shortName = 'f' );
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
