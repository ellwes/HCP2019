

struct Config {
	


} Config;


void init_solver(const double K, const double heat_start, const double cold_start, const double time) {
 /* Init time */	
 /* Create Cart communicator for NxN processes */
 /* Sub div cart communicator to N row communicator - kanske inte behövs  */
 /* Sub div cart communicator to N col communicator - kanske inte behövs  */
 /* Setup sizes of local matrix tiles */
 /* Create subarray datatype for local matrix tile */
 /* Create data array */
 /* Init data array with heat_start and cold_start */

}




void cleanup() {
 /* Set file view and write to file */
 /* Close file */
}





void calc_heat() {
 
 for(double t = config.time; t > 0; t--) {
 	/* Broadcast values to close cells */
 	/* Calculate current theta */
 }

}
