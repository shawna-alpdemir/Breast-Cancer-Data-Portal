function table_to_csv(source) {
    var columns = Object.keys(source.data)
    var nrows = source.get_length()
    var lines = [columns.join(',')]

    for (let i = 0; i < nrows; i++) {
        let row = [];
        for (let j = 0; j < columns.length; j++) {
            var column = columns[j]
            row.push(source.data[column][i].toString())
        }
        lines.push(row.join(','))
    }

    return lines.join('\n').concat('\n')
}

var filename = 'Protein_Summary_Correlation_Table.csv'
var filetext = table_to_csv(source)
var blob = new Blob([filetext], { type: 'text/csv;charset=utf-8;' })

//addresses IE
if (navigator.msSaveBlob) {
    navigator.msSaveBlob(blob, filename)
} else {
    var link = document.createElement('a')
    link.href = URL.createObjectURL(blob)
    link.download = filename
    link.target = '_blank'
    link.style.visibility = 'hidden'
    link.dispatchEvent(new MouseEvent('click'))
}