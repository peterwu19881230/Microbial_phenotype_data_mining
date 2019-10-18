<?php
/*
This code will be used to insert annotations into the website after strain pages are created.

Execute by: php NicholsStrainsCreation.php -w /Library/WebServer/Documents/omp/peter/ -t page_title -f output.csv
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
    
    

			$page_title = $this->getOption('page_title');

            #once I have the paper I am looking for, I search by allele, which in this case is also the strain name
            $dbr = wfGetDB( DB_SLAVE );
            $result = $dbr->select(
                array('page'),
                '*',
                array("page_title LIKE '%$page_title%' "),
                __METHOD__,
                array()
            );

			
			$output_csv = $this->getOption('output_csv');
            $fp = fopen($output_csv, 'w');
                
			#for each strain page found:
            foreach($result as $x){

                #Uses the pages found and creates a title object and finds the tables on the strain pages that were previously created using the pombe_knockout_insert.php code
                $newtitle = Title::newFromText($x->page_title);
                $new_page = new WikiPageTE($newtitle);
                echo $newtitle."\n";



                #Get annotation rows (based on this: https://github.tamu.edu/HuLab/wiki-extensions-mods/blob/79c95b9fdaa7b5f381d97e1fef5c14c2e5411d77/GO_nr/GO_nr.php)
                $newtable = $new_page->getTable("OMP_annotation_table");
                $box = $newtable[0];
                
                
                
                
                
                foreach($box->rows as $row){
                	$data = explode('||', $row->row_data); 
                	fputcsv($fp, $data); 
                }

				


            }

			fclose($fp);

        
    }
    private function parse_parameters(){ #Ref: https://www.mediawiki.org/wiki/Manual:Writing_maintenance_scripts
        $this->addOption( "page_title", "put your page title", $required = true, $withArg = true, $shortName = 't' );
        $this->addOption( "output_csv", "put your output csv name", $required = true, $withArg = true, $shortName = 'f' );
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
