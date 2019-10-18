<?php
/*
Correct strain names by search and replace. (1. title change: JWXXXX-number -> JWXXXX 2. Strain info table: strain name: JWXXXX-number ->JWXXXX, synonyms: JWXXXX-number -> JWXXXX)

(!) After correcting the pages, they no longer appear in the category page (). I think there will be a way to fix that so I still ran this code
(!)This page has been corrected before running the script (added: coli): OMP ST:1173 ! Escherichia K-12 JW5596-1
(!)There were 2 pages for JW5249: JW5249 and JW5249-1. I deleted JW5249-1 and keep JW5249. Also I added JW5249-1 as a synonym in the JW5249 page

Test by: php JWstrainSearchAndReplace.php -w /Library/WebServer/Documents/omp/sandy/

The test website's URL: https://microbialphenotypes.org/sandy/index.php/Main_Page

All the existing JW pages: https://microbialphenotypes.org/sandy/index.php/Category:OMP_ST:800_!_Escherichia_coli_K-12_BW25113_derivatives
*/
$params = getopt( "w:" );
set_IP($params['w']);
$maintClass = "NicholsMaint";
class NicholsMaint extends Maintenance {
    public function __construct() {
        parent::__construct();
        $this->parse_parameters();
    }

    public function execute() {


            #once I have the paper I am looking for, I search by allele, which in this case is also the strain name
            $dbw = wfGetDB( DB_SLAVE );
            $result = $dbw->select(
                array('page'),
                '*',
                array("page_title REGEXP 'Escherichia_coli_K-12_JW[0-9]{4}-'"), #There exists redirects like JW 5335-1 , JW659-5. I don't want to find them.
                __METHOD__,
                array()
            );


            $num=0;
            foreach($result as $x) { #If no pages are found, this foreach() is not run. I don't know why



                #Uses the pages found and creates a title object and finds the tables on the strain pages that were previously created using the pombe_knockout_insert.php code
                $title = Title::newFromText($x->page_title);
                echo $title."\n";
                $old_page = new WikiPageTE($title);
                $oldtable = $old_page->getTable("Strain_info_table");
                if(!isset($oldtable[0])){ echo "->This page has been redirected\n"; continue;} #I don't want the re-directed empty page be processed. Otherwise it will cause Fatal Error


                if(strpos($title, "yraP(del)::FRT-cat-FRT") !== false){ #These are the weird pages I don't understand and want to escape. Ref: https://stackoverflow.com/questions/4366730/how-do-i-check-if-a-string-contains-a-specific-word
                    echo $title,"\n -> is not what I want and will be neglected"."\n";
                    continue;
                }



                echo " ->This page contains dash and number and should be replaced"."\n";
                $num++;

                //Edit and replace the title


                #Check if the page has been redirected (if so the strain info table should be empty)
                $newtitle=preg_replace('/-[0-9]{1}$/',"",$title); echo "Correction: ".$newtitle."\n"; #-[0-9]{1}$ needs to be surrounded by /. Ref:  https://phpenthusiast.com/blog/php-regular-expressions
                $t = Title::newFromText($newtitle);





                #Change the page title
                $context = new RequestContext();
                $movepage = $old_page->move($t, $context);


                //Verify that the updated title is in the database (there will be 2 titles, the additional one which has the old title is there for redirection):

                # login to SQL and search for that page:
                # login on Tetramer
                # mysql -u peterwu => enter password
                # use omp_sandy_wikidb
                # SELECT * FROM page WHERE page_title LIKE '%ST:983_!%';


                $wikiPage = new WikiPageTE($t);  # WikiPageTE is a class defined in WikiPageTE.php
                $newtable = $wikiPage->getTable("Strain_info_table");
                //Remove the -number in the Strain Name box in the Strain Info Table
                $box = $newtable[0];



                # I don't think this is that useful: https://github.tamu.edu/HuLab/wiki-maintenance/blob/master/wiki_editors/SearchReplace.php
                # Will this be the solution?: https://github.tamu.edu/HuLab/wiki-maintenance/blob/master/wiki_editors/EditTablesByTemplate.php


                $strainInfo=$box->rows[0]->row_data_original;


                list($strainName,$synonyms,$TaxonInfo,$Genotype,$StrainRef,$StrainAvail,$ancestry,$annotatedPheno)=explode("||",$strainInfo);
                $newStrainName=preg_replace('/-[0-9]{1}$/',"",$strainName);
                preg_match('/JW[0-9]{4}-[0-9]{1}$/',$title,$matches); #search for a pattern and store it into the variable $matches
                $oldStrainName=$matches[0];
                $newSynonyms=$oldStrainName;



                $box->delete_row(0); #Note: if I delete the row and try to revert it using the webpage, I have to edit the table and save it, otherwise the reverted data will not be put back into the SQL
                #I verified with Sandy that doing this is fine. No need for additional acts to delete things in the SQL
                # 0 inside delete_row is the row index  (in the strain info table the 1st column (not counting the name column) is considered the 1st row)


                $box->insert_row(
                    $newStrainName."||".
                    $newSynonyms."||".
                    $TaxonInfo."||".
                    $Genotype."||".
                    $StrainRef."||".
                    $StrainAvail."||".
                    $ancestry."||".
                    $annotatedPheno
                );



                $wikiPage->touch();








#if($num==1)break;  #This line is used to test the code: $num==1: only 1 page will be corrected
            }

        echo $num." of page(s) has/have been processed"."\n";

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