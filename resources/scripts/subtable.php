<?php
$dataFile = $_POST['dataFile'];
$ids = explode(",",'0,'.$_POST['ids']);
$tableFile = "subtable.xls";

header("Content-type: text/plain");
header("Content-Disposition: attachment; filename=$tableFile");


//unlink($tableFile) or die("Couldn't delete file");
//print_r($ids);
//die;

$fn = fopen($dataFile,"r");
//$fp = fopen($tableFile, 'w') or die('I cant write');

$i = 0;
while(! feof($fn))  {
  $result = fgets($fn);
  
  if (in_array($i,$ids)){
    //fwrite($fp, $result);
    print $result;
  }
  $i++;
}
fclose($fn);
//fclose($fp);
//die;
//header("Location: ".$tableFile);

exit();
?>
