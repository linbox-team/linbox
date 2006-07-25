/* matrix-stream-readers.h
 * Here is where all the formats (each of which is a subclass of
 * MatrixStreamReader) are defined, in two places:
 *
 * First, the macro __MATRIX_STREAM_READERDEFS is put in the init() function of
 * MatrixStream and it creates all the format readers.  For each format reader,
 * a line of this macro should read: addReader( new MyReaderType() );
 *
 * Second, so those statements actually compile, the file containing each format
 * reader should be included with a line of the form: #include "my-reader.h"
 */

#define __MATRIX_STREAM_READERDEFS \
	addReader( new MatrixMarketReader<Field>() );		\
	addReader( new MapleDense1Reader<Field>() );		\
	addReader( new MapleSparse1Reader<Field>() );		\
	addReader( new SparseRowReader<Field>() );		\
	addReader( new SMSReader<Field>() );			\
	addReader( new DenseReader<Field>() );

#include "sms.h"
#include "maple-dense1.h"
#include "maple-sparse1.h"
#include "sparse-row.h"
#include "generic-dense.h"
#include "matrix-market.h"
