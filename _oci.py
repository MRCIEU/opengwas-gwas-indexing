import os

import oci


from dotenv import load_dotenv


load_dotenv()


class OCI:
    def __init__(self):
        self.buckets = {
            'data': os.environ['OCI_BUCKET_DATA'],
            'data-chunks': os.environ['OCI_BUCKET_DATA_CHUNKS'],
            'phewas': os.environ['OCI_BUCKET_PHEWAS'],
            'tophits': os.environ['OCI_BUCKET_TOPHITS'],
        }
        config = {
            "user": os.environ['OCI_USER'],
            "fingerprint": os.environ['OCI_FINGERPRINT'],
            "tenancy": os.environ['OCI_TENANCY'],
            "region": os.environ['OCI_REGION'],
            "key_file": os.environ['OCI_KEY_FILE'],
        }
        oci.config.validate_config(config)
        self.object_storage_client = oci.object_storage.ObjectStorageClient(config)
        self.object_storage_client.base_client.session.adapters['http://'].poolmanager.connection_pool_kw['maxsize'] = 200
        self.object_storage_client.base_client.session.adapters['https://'].poolmanager.connection_pool_kw['maxsize'] = 200

    def object_storage_list(self, bucket_key, prefix):
        results = []
        start = ''

        while True:
            r = self.object_storage_client.list_objects(
                namespace_name=os.environ['OCI_NAMESPACE'],
                bucket_name=self.buckets[bucket_key],
                prefix=prefix,
                start=start,
                fields='size',
            )
            results.extend([o.name for o in r.data.objects])
            if not r.data.next_start_with:
                break
            start = r.data.next_start_with

        return results

    def object_storage_upload(self, bucket_key, object_name, object_body):
        return self.object_storage_client.put_object(
            namespace_name=os.environ['OCI_NAMESPACE'],
            bucket_name=self.buckets[bucket_key],
            object_name=object_name,
            put_object_body=object_body,
        )

    def object_storage_delete(self, bucket_key, object_name):
        return self.object_storage_client.delete_object(
            namespace_name=os.environ['OCI_NAMESPACE'],
            bucket_name=self.buckets[bucket_key],
            object_name=object_name
        )

    def object_storage_delete_by_prefix(self, bucket_key, prefix):
        object_names = self.object_storage_list(bucket_key, prefix)
        for o in object_names:
            self.object_storage_delete(bucket_key, o)
        return object_names

    def object_storage_download(self, bucket_key, object_name):
        return self.object_storage_client.get_object(
            namespace_name=os.environ['OCI_NAMESPACE'],
            bucket_name=self.buckets[bucket_key],
            object_name=object_name
        )
