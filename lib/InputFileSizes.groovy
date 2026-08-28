class InputFileSizes {
    Integer input_size
    nextflow.util.MemoryUnit input_size_bytes

    // Constructors
    InputFileSizes(Collection input_files) {
        // input_size = input_files.collect { f -> f.size() }.sum()
        input_size = input_files.inject(0) { sum, f -> sum + f.size() }
        input_size_bytes = input_size.B
    }

    InputFileSizes(nextflow.processor.TaskPath input_file) {
        input_size = input_file.size()
        input_size_bytes = input_size.B
    }

    InputFileSizes(Integer value) {
        input_size = value
        input_size_bytes = input_size.B
    }

    InputFileSizes(Object nullvalue) {
        assert nullvalue == null : "Invalid argument to InputFileSizes()"
        input_size = 0
        input_size_bytes = input_size.B
    }

    InputFileSizes plus(InputFileSizes other) {
        return new InputFileSizes(input_size + other.input_size)
    }

    long getSizeRawBytes() {
        return input_size
    }

    nextflow.util.MemoryUnit getSizeB() {
        return input_size_bytes
    }

    nextflow.util.MemoryUnit getSizeKB() {
        return input_size_bytes.toKilo().KB
    }

    nextflow.util.MemoryUnit getSizeMB() {
        return input_size_bytes.toMega().MB
    }

    nextflow.util.MemoryUnit getSizeGB() {
        return input_size_bytes.toGiga().GB
    }
}