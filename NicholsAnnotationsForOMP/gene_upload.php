<?php
 #Goal: upload the formatted_unique_gene_names.csv made from Nichols' gene names
 #The Database and the table have to be manually created first.Excecute by: php (this file) -f (your csv)
 #This file was modified from Sandy's
 
  
 $params = getopt( "w:" );
 set_IP($params['w']);
 /*
 Class definition for the desired maintenance object.  
 This must be defined BEFORE the execution lines at the end of the file:
 
 */
 
 $maintClass = "gene_upload";
 
 class gene_upload extends Maintenance {
 
 	public function __construct() {
 		parent::__construct();
 		$this->parse_parameters();
 	}
 
 	/*
 	The guts of the object... do the maintenance task;
 	Replace the content of this function with your code
 	*/
 	public function execute() {
 		# example of using a command line argument
 		$filename = $this->getOption('file');
 		echo "opening $filename\n";
 		$fh = fopen($filename, 'r');


 		for ($i=0; $i<3931;$i++){
 			$line = fgets($fh);
 			#echo "$i $line";
 			
 			if ($line{0} == "\r\n" ||trim($line) == ""){
 				continue;
 			}
 	
 			list ($gene)= explode(",",trim($line));
 			echo "$gene\n";
 			
 			
 			
 			$dbw = wfGetDB( DB_MASTER );
 			$result = $dbw->insert(
 				'Nichols_genes.gene_names',
 				array(
 					'geneID' => null,
 					'gene_name' => $gene,
 				)
 			); 
 			
 			
 		} #close loop 
 	}
   
 	/*
 	
 	to get parameters from the shell add a block of lines of the form
 	$this->addOption( $name, $description, $required = false, $withArg = false, $shortName = false ) 
 	
 	This allows us to use the getter method $this->getOption($name);, and also automatically adds to the help text
 	See maintenance/Maintenance.php for what the arguments mean.
 	
 	*/
 	private function parse_parameters(){
 		$this->addOption( "file", "this is the file we are opening", $required = false, $withArg = true, $shortName = 'f' );	
 	}
 }
 
 require_once( RUN_MAINTENANCE_IF_MAIN );
 
 /*
 Function to set the global variable $IP and include the abstract
 class Maintenance
 */
 function set_IP($path){
    global $IP;
    if ( isset($path) && is_file("$path/maintenance/Maintenance.php") ){
        $IP = $path;
        require_once( $IP . "/maintenance/Maintenance.php" );
        return $path;
    } else {

        echo "Default path to wiki maintenance folder was not specified. Press Enter to continue by using the default path, Enter ant other vavalues to quit:\n";
        echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
        echo "/Library/WebServer/Documents/omp/peter/maintenance/Maintenance.php\n";
        echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

        $handle = fopen ("php://stdin","r");
        $line = fgets($handle);
        if($line=="\n"){
            fclose($handle);
            //I set my wiki as the default path 
            require_once( "/Library/WebServer/Documents/omp/peter/maintenance/Maintenance.php" );
            echo "Default path was used.\n";

        }
        else{
            echo "ABORTING!\n";
            die;
        }

    }

}